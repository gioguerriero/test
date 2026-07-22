function [c, ceq] = nonlinear_constraints_moon2comet_tof(x, c_param, target_pos_synodic, new_tof, theta_fixed, scaling, dv_max)
% nonlinear_constraints_moon2comet_tof  Moon-to-comet constraints with the Moon
% phase angle and time of flight fixed; enforces arrival at the comet position.
%
% Inputs:
%   x                  - design variables [Vinf fpa oop beta dvx dvy dvz]
%   c_param            - constants struct
%   target_pos_synodic - target comet position in the synodic frame [1x3]
%   new_tof            - fixed time of flight [non-dimensional]
%   theta_fixed        - fixed Moon phase angle [rad]
%   scaling            - scaling factor applied to the DSM delta-v variables
%   dv_max             - maximum allowed delta-v (used by the commented constraint)
%
% Outputs:
%   c   - inequality constraints (<= 0)
%   ceq - equality constraints (= 0), comet position match

    % c(x)   <= 0  -> inequalities
    % ceq(x) = 0   -> equalities

    %% Unpack decision variables
    Vinf  = x(1);
    theta = theta_fixed;
    fpa = x(2);
    out_of_plane  = x(3);
    tof   = new_tof;
    beta  = x(4);
    dv = [x(5); x(6); x(7)]./(scaling);

    %% Initialize outputs
    c   = [];
    ceq = [];

    % The only required constraint is arrival at the target position

    %% Integrate from the Moon flyby up to the DSM
    synodic_moon = moon_state(theta, c_param);
    vinf_vec = vinf_rotation(synodic_moon', Vinf, fpa, out_of_plane);
    CI_initial_state = [synodic_moon(1:3)'; synodic_moon(4:6)' + vinf_vec./c_param.Vstar];

    tspan = [0 tof*(beta)]; 

    if tspan(2)>(3600/c_param.Tstar)
        opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
        [t_moon2DSM, S_moon2DSM] = ode45(@(t,S) CR3BP(t,S,c_param.mu), tspan, CI_initial_state, opt);
        pre_DSM_state = S_moon2DSM(end,:)';
        pre_DSM_time = t_moon2DSM(end);
    else
        pre_DSM_state = CI_initial_state;
        pre_DSM_time = 0;
    end

    % Apply the DSM maneuver
    post_DSM_state = [pre_DSM_state(1:3); pre_DSM_state(4:6) + dv];

    % Re-propagate to the comet
    tspan = [0 tof*(1-beta)];

    if tspan(2)>(3600/c_param.Tstar)
        opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
        [t_DSM2comet, S_DSM2comet] = ode45(@(t,S) CR3BP(t,S,c_param.mu), tspan, post_DSM_state, opt);
        arrival_state = S_DSM2comet(end,:)';
    else
        arrival_state = post_DSM_state;
    end

    ceq = [ceq; (arrival_state(1:3) - target_pos_synodic).*1e+2];
    % c = [c; norm(dv) - dv_max]; % dv in questo caso è adimensionale, non metto la cost function ma metto un vincolo che deltaV non deve essere troppo più grosso della guess (di 20m/s?)

end