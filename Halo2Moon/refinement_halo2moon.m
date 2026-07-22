function [refined_traj, tof_halo2moon, exitflag] = refinement_halo2moon(c, solution, dv_max, vinf_fmincon, fpa_fmincon, out_of_plane_fmincon, h_min_earth, eps_vel_ms, S_halo, unstable_dir, fmincon_vinf_results, min_days_between)
% refinement_halo2moon  fmincon refinement of a grid-search Halo-to-Moon
% transfer, minimizing the total time of flight under closure constraints.
%
% Inputs:
%   c                    - constants struct
%   solution             - grid-search match struct (angles, TOFs, moon_synodic)
%   dv_max               - maximum allowed delta-v [km/s]
%   vinf_fmincon         - reference v-infinity magnitude [km/s]
%   fpa_fmincon          - reference flight-path angle [rad]
%   out_of_plane_fmincon - reference out-of-plane angle [rad]
%   h_min_earth          - minimum Earth altitude constraint [km]
%   eps_vel_ms           - injection velocity perturbation [m/s]
%   S_halo               - halo orbit states
%   unstable_dir         - unstable eigenvectors
%   fmincon_vinf_results - reference v-infinity results struct
%   min_days_between     - min event spacing [days]; only elements 1-2 used here
%
% Outputs:
%   refined_traj  - optimized design variables [fpa oop halo_theta tof_moon tof_halo]
%   tof_halo2moon - total Halo-to-Moon time of flight
%   exitflag      - fmincon exit flag

% Set the initial guess from the grid-search solution
moon_fpa = solution.moon_fpa;
moon_out_of_plane = solution.moon_out_of_plane;
halo_theta = solution.halo_theta;
tof_moon = -solution.moon_time_of_encounter;
tof_halo = solution.halo_time_of_encounter;

% If it starts at 0 the design variable never moves; enforce a small minimum
tof_halo = max(tof_halo, 1/c.Tstar);

% Reference V-infinity vector
moon_synodic = solution.moon_synodic; 

v0 = vinf_rotation(moon_synodic, vinf_fmincon, fpa_fmincon, out_of_plane_fmincon);
v0_norm = norm(v0);

% Maximum deflection angle (lunar gravity assist)
muMoon = 4902.800066;   % km^3/s^2
rMoon  = 1737;           % km
hp     = 750;            % km
rM     = rMoon + hp;

delta_max = 2*asin( muMoon ./ (rM * v0_norm^2 + muMoon) );

% Initial guess fmincon
x0 = [moon_fpa, moon_out_of_plane, halo_theta, tof_moon, tof_halo];

% Minimum TOF spacing so the maneuver is not placed right at the Moon flyby or
% at the Halo departure. x(4)=tof_moon (DSM1->flyby), x(5)=tof_halo (inj->DSM1).
lb = [-pi  -pi  -pi   min_days_between(2)*24*3600/c.Tstar   min_days_between(1)*24*3600/c.Tstar];
ub = [pi    pi   pi  1*365*24*3600/c.Tstar  1*365*24*3600/c.Tstar];

% Linear constraints (none)
A = [];
b = [];
Aeq = [];
beq = [];

% Nonlinear constraints (if needed)
nonlcon = @(x) NC_refinement_halo2moon(x, c, dv_max, h_min_earth, delta_max, moon_synodic, v0_norm, eps_vel_ms, S_halo, unstable_dir, fmincon_vinf_results);

% fmincon Options
options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...  % SQP is robust  
    'MaxIterations',70, ...
    'Display','off', ...
    'MaxFunctionEvaluations',400, ... 
    'OptimalityTolerance',1e-5, ...
    'StepTolerance',1e-10, ...
    'ConstraintTolerance',1e-5);
% 'Display','off', ...
 

% Call fmincon
[x_opt, fval, exitflag, output] = fmincon( ...
    @(x) OF_refinement_moon2comet(x, c), ...
    x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

refined_traj = x_opt;
tof_halo2moon = fval;



end