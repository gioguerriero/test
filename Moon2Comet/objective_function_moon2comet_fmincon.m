function dv = objective_function_moon2comet_fmincon(x, c, theta_inertial_final, q, scaling)
% objective_function_moon2comet_fmincon  fmincon cost for the Moon-to-comet
% leg: DSM delta-v (dimensional inertial) plus a time-of-flight penalty.
%
% Inputs:
%   x                    - design variables [Vinf theta fpa oop tof beta dvx dvy dvz]
%   c                    - constants struct
%   theta_inertial_final - synodic frame angle at the comet epoch [rad]
%   q                    - weight on time of flight in the cost
%   scaling              - scaling factor applied to the DSM delta-v variables
%
% Outputs:
%   dv - scalar cost (DSM delta-v [km/s] + q*tof)

    %% Compute the delta-v in the inertial frame (minimized in dimensional inertial units)

    Vinf  = x(1);
    theta = x(2);
    fpa = x(3);
    out_of_plane  = x(4);
    tof   = x(5);
    beta  = x(6);

    synodic_moon = moon_state(theta, c);
    vinf_vec = vinf_rotation(synodic_moon', Vinf, fpa, out_of_plane);
    CI_initial_state = [synodic_moon(1:3)'; synodic_moon(4:6)' + vinf_vec./c.Vstar];

    tspan = [0 tof*(beta)]; 

    if tspan(2)>(3600/c.Tstar)
        opt = odeset('AbsTol',1e-8,'RelTol',1e-8);
        [t_moon2DSM, S_moon2DSM] = ode45(@(t,S) CR3BP(t,S,c.mu), tspan, CI_initial_state, opt);
        pre_DSM_state = S_moon2DSM(end,:)';
        pre_DSM_time = t_moon2DSM(end);
    else
        pre_DSM_state = CI_initial_state;
        pre_DSM_time = 0;
    end

    post_DSM_state = [pre_DSM_state(1:3); pre_DSM_state(4:6) + [x(7); x(8); x(9)]./(scaling)];

    theta_initial = theta_inertial_final - tof;
    pre_DSM_state_cartesian = synodic2car(pre_DSM_state', pre_DSM_time, c.mu, theta_initial);
    post_DSM_state_cartesian = synodic2car(post_DSM_state', pre_DSM_time, c.mu, theta_initial);

    dv = norm(post_DSM_state_cartesian(4:6) - pre_DSM_state_cartesian(4:6))*c.Vstar + q*tof;


end
