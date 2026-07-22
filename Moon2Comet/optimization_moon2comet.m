function [result_ga, result_fmincon, result_fmincon_modtof] = optimization_moon2comet(c, target_pos_inertial, vinf_data, theta_inertial_final, target_pos_synodic, epoch_encounter, q, min_days_between)
% optimization_moon2comet  Three-stage optimizer for the Moon-to-comet leg:
% GA global search, fmincon refinement, then a fixed-TOF fmincon pass.
%
% Inputs:
%   c                    - constants struct
%   target_pos_inertial  - comet position, inertial [km]
%   vinf_data            - struct with vinf lower_bound/upper_bound [km/s]
%   theta_inertial_final - synodic frame angle at the comet epoch [rad]
%   target_pos_synodic   - comet position in the synodic frame [1x3]
%   epoch_encounter      - comet encounter epoch [ET s]
%   q                    - weight on time of flight in the cost
%   min_days_between     - min spacing [inj->DSM1, DSM1->flyby, flyby->DSM2, DSM2->comet] [days]
%
% Outputs:
%   result_ga             - GA result (x_opt, fval)
%   result_fmincon        - fmincon output struct (stage 2)
%   result_fmincon_modtof - fixed-TOF fmincon output struct (stage 3)
%
% Uses spacing elements 3 (flyby->DSM2 = beta*tof) and 4 (DSM2->comet = (1-beta)*tof):
%   - stage 2 (fmincon, variable tof): nonlinear constraint in nonlinear_constraints_moon2comet
%   - stage 3 (fmincon_tof, fixed tof): bound on beta
D_flyby_dsm2_adim = min_days_between(3) * 24*3600 / c.Tstar;
D_dsm2_comet_adim = min_days_between(4) * 24*3600 / c.Tstar;

% Start from a lunar orbit and optimize, in CR3BP dynamics, the parameters
% that minimize the DSM delta-v:
%   - Vinf
%   - crank angle
%   - pump angle
%   - total time of flight
%   - beta (beta * tof defines the epoch at which the DSM is performed)



%% 1) genetic algorithm (ga)

    % Number of decision variables
    nvars = 6;

    % Bounds (note: the beta lower bound should avoid a DSM right after the flyby)
    % [Vinf, position in moon orbit, flight path angle, out of plane angle, TOF, beta]
    % zero of "position in moon orbit" points toward Sun-Earth L2, increasing counter-clockwise
    lb = [ vinf_data.lower_bound, -pi   -pi,    -pi,   1e-4,    1e-2 ];
    ub = [ vinf_data.upper_bound,  pi,    pi,     pi,  4*365*24*3600/c.Tstar,    1.0-1e-2 ];

    % Linear constraints (none for now)
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    % Nonlinear constraints (if needed)
    nonlcon = [];

    % GA options
    options = optimoptions('ga',...
        'PopulationSize', 1000,...
        'MaxGenerations', 100,...
        'Display','iter',...
        'UseParallel', true,...      
        'FunctionTolerance', 1e-5);

    % Call GA
    [x_opt, fval] = ga(@(x) objective_function_moon2comet_ga(x, c, target_pos_inertial, theta_inertial_final, q), ...
                       nvars, A, b, Aeq, beq, lb, ub, nonlcon, options);

    % Output
    result_ga.x_opt = x_opt;
    result_ga.fval  = fval;


%% 2) gradient based refinement algorithm (fmincon) + 2BP -> CR3BP dynamics shift (for post DSM)

% Design variables: the previous ones plus the DSM direction (2 angles) and magnitude

% Bounds
lb = [ vinf_data.lower_bound,  -pi,  -pi,  -pi,  1e-4,                  1e-2,  -inf, -inf, -inf];
ub = [ vinf_data.upper_bound,   pi,   pi,   pi,  4*365*24*3600/c.Tstar,  1.0-1e-2,   inf,  inf,  inf];

% Linear constraints (none)
A = [];
b = [];
Aeq = [];
beq = [];

DSM_info = get_DSM_info(result_ga.x_opt, c, theta_inertial_final, target_pos_inertial);
% Scale DSM_info to order 1: normalise it so its magnitude becomes ~10
scaling1 = 10/norm(DSM_info);

% Initial guess from the genetic algorithm result
x0 = [result_ga.x_opt, DSM_info.*scaling1];

% Nonlinear constraints (if needed)
nonlcon = @(x) nonlinear_constraints_moon2comet(x, c, target_pos_synodic, scaling1, D_flyby_dsm2_adim, D_dsm2_comet_adim);

% fmincon Options
options = optimoptions('fmincon', ...
    'Display','iter', ...
    'Algorithm','sqp', ...              % SQP is robust
    'MaxIterations',500, ...
    'MaxFunctionEvaluations',5000, ...
    'OptimalityTolerance',1e-6, ...
    'StepTolerance',1e-10, ...
    'ConstraintTolerance',1e-7);


% Call fmincon
[x_opt, fval, exitflag, output] = fmincon( ...
    @(x) objective_function_moon2comet_fmincon(x, c, theta_inertial_final, q, scaling1), ...
    x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

% Output
result_fmincon = output;


%% Time of flight adjustment
% Same problem as stage 2, but now the TOF becomes a constant. It is adjusted
% based on the actual Moon position at the flyby (increased or decreased as needed).

% Compute the Moon position at the flyby epoch to correct the TOF
moon_flyby_epoch = epoch_encounter - x_opt(5)*c.Tstar;

[state_earth, ~] = cspice_spkezr('EARTH', moon_flyby_epoch, 'ECLIPJ2000', 'NONE', 'SUN');
[state_moon, ~] = cspice_spkezr('MOON', moon_flyby_epoch, 'ECLIPJ2000', 'NONE', 'SUN');

pos_earth = state_earth(1:3);
pos_moon = state_moon(1:3);

% Moon angle w.r.t. Earth, measured from the Sun->Earth axis (x-y plane only)
r_moon_from_earth = pos_moon(1:2) - pos_earth(1:2);   % Earth->Moon vector (x,y)

% Reference unit vector: Sun->Earth direction
r_sun_earth_hat      = pos_earth(1:2) / norm(pos_earth(1:2));
% Perpendicular unit vector (90 deg counter-clockwise) to build the rotated frame
r_sun_earth_hat_perp = [-r_sun_earth_hat(2); r_sun_earth_hat(1)];

% Projection of Earth->Moon onto the (ref, perp) frame
x_proj = dot(r_moon_from_earth, r_sun_earth_hat);
y_proj = dot(r_moon_from_earth, r_sun_earth_hat_perp);

% Angle in [-pi, pi]: zero on the Sun->Earth axis, positive counter-clockwise
theta_moon_real = atan2(y_proj, x_proj);   % [rad]
% Wrap the angle into [0, 2pi]
if theta_moon_real<0
    theta_moon_real = theta_moon_real + 2*pi;
end

% Difference with the computed encounter angle
theta_moon_computed = x_opt(2); % also in [-pi, pi]
if theta_moon_computed<0
    theta_moon_computed = theta_moon_computed + 2*pi;
end

% Raw difference in (-2pi, 2pi): positive = Moon ahead, negative = Moon behind
relative_angle_raw = theta_moon_real - theta_moon_computed;

relative_angle = mod(relative_angle_raw + pi, 2*pi) - pi;  % [rad]

% n_moon = 2.6617e-06; % moto medio (siderale) della Luna intorno alla Terra [rad/s]
n_moon = 2*pi/(29.530589*86400);   % SYNODIC mean motion [rad/s]
dt = relative_angle / n_moon;  % time correction [s]

% Increase or decrease the TOF while keeping the Moon encounter angle fixed
new_tof = x_opt(5) + dt/c.Tstar;
theta_fixed = x_opt(2);


% Bounds -> the TOF is no longer a design variable.
% M2C spacing enforced as a bound on beta: new_tof is fixed, so
%   beta*new_tof >= D_flyby_dsm2  ->  beta >= D_flyby_dsm2/new_tof
%   (1-beta)*new_tof >= D_dsm2_comet  ->  beta <= 1 - D_dsm2_comet/new_tof
beta_lb = max(1e-2,     D_flyby_dsm2_adim / new_tof);
beta_ub = min(1.0-1e-2, 1 - D_dsm2_comet_adim / new_tof);
if beta_lb > beta_ub
    warning(['optimization_moon2comet: M2C spacing infeasible for TOF=%.1f d ' ...
             '(flyby->DSM2 %.1f + DSM2->comet %.1f days). Relaxing to the original bounds.'], ...
            new_tof*c.Tstar/86400, min_days_between(3), min_days_between(4));
    beta_lb = 1e-2;  beta_ub = 1.0-1e-2;
end
lb = [ vinf_data.lower_bound,  -pi,  -pi,  beta_lb,  -inf, -inf, -inf];
ub = [ vinf_data.upper_bound,   pi,   pi,  beta_ub,   inf,  inf,  inf];

% Linear constraints (none)
A = [];
b = [];
Aeq = [];
beq = [];

deltaV_adim = x_opt(7:9)./scaling1;
scaling2 = 1/(norm(deltaV_adim));

% Initial guess from the fmincon result, dropping the tof
x0 = [x_opt(1), x_opt(3:4), x_opt(6), deltaV_adim.*scaling2];

% Add a constraint on dv_max (no cost function, since very small delta-v tends not to converge) -> currently inactive
dv_margin = 30 / (c.Vstar*1e3); % [m/s]
dv_max = norm(deltaV_adim) + dv_margin;

% Nonlinear constraints (if needed)
nonlcon = @(x) nonlinear_constraints_moon2comet_tof(x, c, target_pos_synodic, new_tof, theta_fixed, scaling2, dv_max);

% fmincon Options
options = optimoptions('fmincon', ...
    'Display','iter', ...
    'Algorithm','sqp', ...              % SQP is robust
    'MaxIterations',500, ...
    'MaxFunctionEvaluations',5000, ...
    'OptimalityTolerance',1e-6, ...
    'StepTolerance',1e-10, ...
    'ConstraintTolerance',1e-5);   % stage 3: ~15 km on the comet match (CR3BP, later refined in ephemeris)

% Call fmincon
[x_opt, fval, exitflag, output] = fmincon( ...
    @(x) objective_function_moon2comet_tof(x, c, theta_inertial_final, new_tof, theta_fixed, scaling2), ...
    x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

% Output
output.bestfeasible.x = [output.bestfeasible.x(1), theta_fixed, output.bestfeasible.x(2:3), new_tof, output.bestfeasible.x(4), output.bestfeasible.x(5:7)./(scaling2)];
result_fmincon_modtof = output;




end

