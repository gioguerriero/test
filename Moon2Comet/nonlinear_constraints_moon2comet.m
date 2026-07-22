function [c, ceq] = nonlinear_constraints_moon2comet(x, c_param, target_pos_synodic, scaling, D_flyby_dsm2_adim, D_dsm2_comet_adim)
% nonlinear_constraints_moon2comet  Constraints for the Moon-to-comet leg:
% arrival at the target comet position and minimum maneuver spacing.
%
% Inputs:
%   x                  - design variables [Vinf theta fpa oop tof beta dvx dvy dvz]
%   c_param            - constants struct
%   target_pos_synodic - target comet position in the synodic frame [1x3]
%   scaling            - scaling factor applied to the DSM delta-v variables
%   D_flyby_dsm2_adim  - minimum flyby->DSM2 time [non-dimensional]
%   D_dsm2_comet_adim  - minimum DSM2->comet time [non-dimensional]
%
% Outputs:
%   c   - inequality constraints (<= 0)
%   ceq - equality constraints (= 0), comet position match

    % c(x)   <= 0  -> inequalities
    % ceq(x) = 0   -> equalities

    %% Unpack decision variables
    Vinf  = x(1);
    theta = x(2);
    fpa = x(3);
    out_of_plane  = x(4);
    tof   = x(5);
    beta  = x(6);
    dv = [x(7); x(8); x(9)];

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
    post_DSM_state = [pre_DSM_state(1:3); pre_DSM_state(4:6) + [x(7); x(8); x(9)]./(scaling)];

    % Re-propagate to the comet
    tspan = [0 tof*(1-beta)];

    if tspan(2)>(3600/c_param.Tstar)
        opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
        [t_DSM2comet, S_DSM2comet] = ode45(@(t,S) CR3BP(t,S,c_param.mu), tspan, post_DSM_state, opt);
        arrival_state = S_DSM2comet(end,:)';
    else
        arrival_state = post_DSM_state;
    end

    ceq = [ceq; arrival_state(1:3) - target_pos_synodic];

    %% M2C maneuver spacing (tof = x(5) is variable -> nonlinear constraint)
    %   flyby->DSM2 = beta*tof       >= D_flyby_dsm2_adim
    %   DSM2->comet = (1-beta)*tof   >= D_dsm2_comet_adim
    c = [c;
         D_flyby_dsm2_adim - tof*beta;
         D_dsm2_comet_adim - tof*(1-beta)];

end