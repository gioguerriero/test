function manifold = build_halo_manifold_db(S_halo, unstable_dir, c, N, eps_vel_ms, tmax_years, h_min_earth, points)
% BUILD_HALO_MANIFOLD_DB  Propagate unstable manifold from N points on the halo orbit.
%
%   Inputs:
%       S_halo        - halo orbit states [Npts x 6]
%       unstable_dir  - unstable directions [Npts x 3]
%       c             - constants structure
%       N             - number of equally spaced departure points
%       eps_vel_ms    - velocity perturbation amplitude [m/s]
%       tmax_years    - propagation time [years]
%       h_min_earth   - minimum altitude above Earth [km]
%       points        - number of samples per propagated arc
%
%   Output:
%       manifold      - struct array with only valid trajectories

if nargin < 4 || isempty(N),            N = 100;            end
if nargin < 5 || isempty(eps_vel_ms),   eps_vel_ms = 10;    end
if nargin < 6 || isempty(tmax_years),   tmax_years = 1;     end
if nargin < 7 || isempty(h_min_earth),  h_min_earth = 0;    end

% Earth position in CR3BP (normalized)
r_earth = [1 - c.mu, 0, 0];

% Convert altitude threshold from km to DU
h_min_DU = h_min_earth / c.Lstar;

% Earth's radius in DU
r_earth_radius_DU = c.rEarth / c.Lstar;

r_min_allowed = r_earth_radius_DU + h_min_DU;

% Equally spaced indices on halo
idx = round(linspace(1, size(S_halo,1), N));

% Propagation time span
t_end = tmax_years*365*24*3600/c.Tstar;
tspan = linspace(0, t_end, points);

% Maximum distance from Earth: 5 Earth-Moon distances
r_max = 10 * c.rMoon_ad;  % in DU
opt = odeset('AbsTol',1e-8, 'RelTol',1e-8, ...
             'Events', @(t,S) escape_event(t, S, r_earth, r_max));

% Velocity perturbation (adimensional)
eps_vel = eps_vel_ms / (1e3 * c.Vstar);

manifold = struct([]);
count = 0;

for k = 1:N
    k
    i = idx(k);

    x0 = S_halo(i,:)';

    vu = unstable_dir(i,:)';
    vu = vu / norm(vu);

    x0_pert = [x0(1:3);
               x0(4:6) + eps_vel*vu];

    [t_vec, S_vec] = ode45(@(t,S) CR3BP(t,S,c.mu), ...
                       tspan, x0_pert, opt);

    % ======================================================
    % CHECK MIN DISTANCE FROM EARTH
    % ======================================================

    r_vec = S_vec(:,1:3);

    % distance from the Earth centre along the whole trajectory
    dist_earth = vecnorm(r_vec - r_earth, 2, 2);

    % ======================================================
    % REFINED MINIMUM DISTANCE CHECK (local refinement)
    % ======================================================

    % Find index of minimum distance in the coarse trajectory
    [~, idx_min] = min(dist_earth);
    S_min = S_vec(idx_min,1:6);

    % Define ±5 days window (in nondimensional time)
    delta_t_days = 5;
    delta_t_nd = (delta_t_days * 24 * 3600) / c.Tstar;

    % Refined time grid (1000 points)
    t_refined_post = linspace(0, delta_t_nd, 500);
    t_refined_pre = linspace(0, -delta_t_nd, 500);

    % Re-integrate trajectory locally with higher resolution
    [~, S_refined_post] = ode45(@(t,S) CR3BP(t,S,c.mu), ...
                           t_refined_post, S_min, opt);

    [~, S_refined_pre] = ode45(@(t,S) CR3BP(t,S,c.mu), ...
                           t_refined_pre, S_min, opt);

    % Recompute distances from Earth on refined trajectory
    r_refined_post = S_refined_post(:,1:3);
    r_refined_pre = S_refined_pre(:,1:3);
    dist_refined_post = vecnorm(r_refined_post - r_earth, 2, 2);
    dist_refined_pre = vecnorm(r_refined_pre - r_earth, 2, 2);
    dist_refined = [dist_refined_pre; dist_refined_post];

    % Updated minimum distance
    d_min = min(dist_refined);

    % Filter
    if d_min < r_min_allowed
        continue; % discard trajectory
    end

    % passed the filter -> save
    count = count + 1;

    manifold(count).theta   = atan2(x0(3), x0(2));
    manifold(count).x0      = x0;
    manifold(count).x0_pert = x0_pert;
    manifold(count).t       = t_vec';
    manifold(count).state   = S_vec;
    manifold(count).dmin    = d_min;
end

% trim any empty cells
manifold = manifold(1:count);

disp("Manifold propagation completed");
fprintf("Kept %d / %d trajectories after Earth proximity filter\n", count, N);

save("manifold_db.mat", "manifold")

end


function [value, isterminal, direction] = escape_event(~, S, r_earth, r_max)
% escape_event  ode45 event that stops propagation when the trajectory moves
% beyond r_max from Earth.
    value      = norm(S(1:3)' - r_earth) - r_max;
    isterminal = 1;
    direction  = 1;  % trigger only when moving away
end
