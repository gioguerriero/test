function manifold_back = build_moon_manifold_db(c, theta_moon, vinf_fmincon, ...
    fpa_fmincon, out_of_plane_fmincon, M, tmax_years, h_min_earth, points)
% build_moon_manifold_db  Back-propagate a fan of trajectories departing the
% Moon over a deflection-angle grid, keeping only those clearing a minimum Earth altitude.
%
% Inputs:
%   c                    - constants struct
%   theta_moon           - Moon phase angle in the synodic frame [rad]
%   vinf_fmincon         - reference v-infinity magnitude [km/s]
%   fpa_fmincon          - reference flight-path angle [rad]
%   out_of_plane_fmincon - reference out-of-plane angle [rad]
%   M                    - approximate number of grid points
%   tmax_years           - back-propagation time [years]
%   h_min_earth          - minimum Earth altitude [km]
%   points               - samples per propagated arc
%
% Outputs:
%   manifold_back - struct array of valid back-propagated trajectories

% Defaults
if nargin < 6 || isempty(M),            M = 100;    end
if nargin < 7 || isempty(tmax_years),   tmax_years = 1; end
if nargin < 8 || isempty(h_min_earth),  h_min_earth = 0; end

% Moon state in synodic frame
moon_synodic = moon_state(theta_moon, c)';

% Earth position in CR3BP
r_earth = [1 - c.mu, 0, 0];

% Convert altitude threshold to DU
h_min_DU = h_min_earth / c.Lstar;

% Earth radius in DU
r_earth_radius_DU = c.rEarth / c.Lstar;

r_min_allowed = r_earth_radius_DU + h_min_DU;

% Reference V-infinity vector
v0 = vinf_rotation(moon_synodic, vinf_fmincon, fpa_fmincon, out_of_plane_fmincon);
v0_norm = norm(v0);

% Maximum deflection angle (lunar gravity assist)
muMoon = 4902.800066;   % km^3/s^2
rMoon  = 1737;           % km
hp     = 750;            % km
rM     = rMoon + hp;

delta_max = 2*asin( muMoon ./ (rM * v0_norm^2 + muMoon) );

% Polar uniform grid: equispaced radii, variable angular points per ring
n_r = round(sqrt(M));
r_vec_grid = linspace(0, delta_max, n_r);

% Distribute ~M points: 1 at center, rest proportional to ring index
ring_weights = (1:n_r-1);
pts_per_ring = round((M - 1) * ring_weights / sum(ring_weights));
pts_per_ring = max(pts_per_ring, 3);

FPA = 0;   % center point
OUT = 0;

for kr = 2:n_r
    n_ang_k = pts_per_ring(kr-1);
    ang_k   = linspace(0, 2*pi, n_ang_k+1); ang_k(end) = [];
    FPA = [FPA; r_vec_grid(kr) * cos(ang_k')];
    OUT = [OUT; r_vec_grid(kr) * sin(ang_k')];
end

M_real = length(FPA);

% Propagation settings (backward in time)
t_end = -tmax_years*365*24*3600/c.Tstar;
tspan = linspace(0, t_end, points);

% Maximum distance from Earth: 5 Earth-Moon distances
r_max = 10 * c.rMoon_ad;  % in DU
opt = odeset('AbsTol',1e-8, 'RelTol',1e-8, ...
             'Events', @(t,S) escape_event(t, S, r_earth, r_max));

manifold_back = struct([]);
count = 0;

for k = 1:M_real
    fpa_k = FPA(k);
    out_k = OUT(k);

    v_inf_k = vinf_rotation(moon_synodic, v0_norm, fpa_k, out_k);

    x0 = moon_synodic;
    x0(4:6) = x0(4:6) + v_inf_k./c.Vstar;

    [t_vec, S_vec] = ode45(@(t,S) CR3BP(t,S,c.mu), ...
                      tspan, x0, opt);

    % ======================================================
    % CHECK MIN DISTANCE FROM EARTH
    % ======================================================

    r_vec = S_vec(:,1:3);

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

    if d_min < r_min_allowed
        continue; % discard trajectory
    end

    % save only if valid
    count = count + 1;

    manifold_back(count).fpa          = fpa_k;
    manifold_back(count).out_of_plane = out_k;
    manifold_back(count).v_inf        = v_inf_k;
    manifold_back(count).x0           = x0;
    manifold_back(count).t            = t_vec';
    manifold_back(count).state        = S_vec;
    manifold_back(count).dmin         = d_min;
    manifold_back(count).moon_synodic = moon_synodic;
end

% Trim struct array
manifold_back = manifold_back(1:count);

disp("Trajectory from moon propagation completed");
fprintf("Kept %d / %d trajectories after Earth proximity filter\n", count, M_real);

end


function [value, isterminal, direction] = escape_event(~, S, r_earth, r_max)
% escape_event  ode45 event that stops propagation when the trajectory moves
% beyond r_max from Earth.
    value      = norm(S(1:3)' - r_earth) - r_max;
    isterminal = 1;
    direction  = 1;  % trigger only when moving away
end
