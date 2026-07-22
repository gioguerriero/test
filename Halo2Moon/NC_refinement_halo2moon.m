function [c, ceq] = NC_refinement_halo2moon(x, c_const, dv_max, h_min_earth, delta_max, moon_synodic, v0_norm, eps_vel_ms, S_halo, unstable_dir, fmincon_vinf_results)
% NC_refinement_halo2moon  Nonlinear constraints for the Halo-to-Moon refinement:
% arc position match, delta-v budget, deflection-angle and Earth-altitude limits.
%
% Inputs:
%   x                    - design variables [fpa oop halo_theta tof_moon tof_halo]
%   c_const              - constants struct
%   dv_max               - maximum allowed delta-v [km/s]
%   h_min_earth          - minimum Earth flyby altitude [km]
%   delta_max            - maximum lunar deflection angle [rad]
%   moon_synodic         - Moon synodic state [1x6]
%   v0_norm              - reference v-infinity magnitude
%   eps_vel_ms           - injection velocity perturbation [m/s]
%   S_halo               - halo orbit states
%   unstable_dir         - unstable eigenvectors
%   fmincon_vinf_results - reference v-infinity results struct
%
% Outputs:
%   c   - inequality constraints (<= 0)
%   ceq - equality constraints (= 0), arc position match

moon_fpa = x(1);
moon_out_of_plane = x(2);
halo_theta = x(3);
tof_moon = x(4);
tof_halo = x(5);

c = [];
ceq = [];

opt = odeset('AbsTol',1e-8,'RelTol',1e-8);

% The constraints are:
% - the two arcs must match in position
% - the departure angle must be within delta_max
% - the delta-v must be below the maximum
% - the closest approach to Earth must be above h_min_earth

%% Integrate the Moon arc
v_inf = vinf_rotation(moon_synodic, v0_norm, moon_fpa, moon_out_of_plane);
S0 = moon_synodic;
S0(4:6) = S0(4:6) + (v_inf)./c_const.Vstar;
tspan = [0 -tof_moon];

[~, S_moon_mat] = ode45(@(t,S) CR3BP(t,S,c_const.mu), tspan, S0, opt);


%% Integrate the Halo arc
[initial_state, idx] = state_finder(halo_theta, S_halo);
eps_vel = eps_vel_ms / (1e3 * c_const.Vstar);
vu = unstable_dir(idx,:)';
vu = vu / norm(vu);

S0_pert = [initial_state(1:3); initial_state(4:6) + eps_vel*vu];

tspan = [0 tof_halo];
[~, S_halo_mat] = ode45(@(t,S) CR3BP(t,S,c_const.mu), tspan, S0_pert, opt);


%% Constraint definition

final_state_moon = S_moon_mat(end,:);
final_state_halo = S_halo_mat(end,:);

% constraints 1 and 2 -> position must match and delta-v must be below dv_max
ceq = [ceq; final_state_halo(1:3)' - final_state_moon(1:3)'];
c = [c; norm(final_state_halo(4:6)-final_state_moon(4:6)) - dv_max/c_const.Vstar];

% The delta angle is obtained from the reference data in fmincon_vinf_results

% constraint 3 -> maximum deflection angle
v_inf_fmincon = vinf_rotation(moon_synodic, fmincon_vinf_results.vinf_fmincon, fmincon_vinf_results.fpa_fmincon, fmincon_vinf_results.out_of_plane_fmincon);
cos_theta = dot(v_inf, v_inf_fmincon) / (norm(v_inf)*norm(v_inf_fmincon));
c = [c; cos(delta_max) - cos_theta];


% constraint 4 -> minimum distance from Earth
% ======================================================
% REFINED MINIMUM DISTANCE CHECK (local refinement) -> for the moment is
% commented (lines 101-123) as it saves computational time
% ======================================================

opt = odeset('AbsTol',1e-8,'RelTol',1e-8);

r_earth = [1 - c_const.mu, 0, 0];

% Combine both arcs into a single state matrix
% Moon arc is back-propagated (times go negative), flip it so time increases
S_full = [flipud(S_moon_mat); S_halo_mat];

% Coarse distance from Earth
dist_earth = vecnorm(S_full(:,1:3) - r_earth, 2, 2);
d_min_coarse = min(dist_earth);

% Skip refinement if coarse minimum is already far from Earth
skip_threshold = (100000 + c_const.rEarth) / c_const.Lstar;

if d_min_coarse > skip_threshold
    % Far enough, use coarse value directly (no extra propagation needed)
    d_min = d_min_coarse;
% else
%     % Close to Earth: refine locally
%     [~, idx_min] = min(dist_earth);
%     S_min = S_full(idx_min, 1:6);
% 
%     % Define +/-5 days window (in nondimensional time)
%     delta_t_nd = (5 * 24 * 3600) / c_const.Tstar;
% 
%     % Refined propagation: forward and backward from closest point
%     t_refined_post = linspace(0, delta_t_nd, 200);
%     t_refined_pre  = linspace(0, -delta_t_nd, 200);
% 
%     [~, S_refined_post] = ode45(@(t,S) CR3BP(t,S,c_const.mu), ...
%                            t_refined_post, S_min, opt);
% 
%     [~, S_refined_pre] = ode45(@(t,S) CR3BP(t,S,c_const.mu), ...
%                            t_refined_pre, S_min, opt);
% 
%     % Minimum distance on refined grid
%     dist_refined = [vecnorm(S_refined_pre(:,1:3)  - r_earth, 2, 2);
%                     vecnorm(S_refined_post(:,1:3) - r_earth, 2, 2)];
% 
%     d_min = min(dist_refined);
end

% Constraint: d_min >= h_min_earth (already adimensional)
c = [c; (h_min_earth + c_const.rEarth)/c_const.Lstar - d_min];


end