function dv = objective_function_moon2comet_tof(x, c, theta_inertial_final, new_tof, theta_fixed, scaling)
% objective_function_moon2comet_tof  Moon-to-comet cost variant with the Moon
% phase angle and time of flight fixed; minimizes the DSM delta-v only.
%
% Inputs:
%   x                    - design variables [Vinf fpa oop beta dvx dvy dvz]
%   c                    - constants struct
%   theta_inertial_final - synodic frame angle at the comet epoch [rad]
%   new_tof              - fixed time of flight [non-dimensional]
%   theta_fixed          - fixed Moon phase angle [rad]
%   scaling              - scaling factor applied to the DSM delta-v variables
%
% Outputs:
%   dv - scalar cost (DSM delta-v [km/s])

    %% Compute the delta-v in the inertial frame (minimized in dimensional inertial units)

    Vinf  = x(1);
    theta = theta_fixed;
    fpa = x(2);
    out_of_plane  = x(3);
    tof   = new_tof;
    beta  = x(4);

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

    post_DSM_state = [pre_DSM_state(1:3); pre_DSM_state(4:6) + [x(5); x(6); x(7)]./scaling];

    theta_initial = theta_inertial_final - tof;
    pre_DSM_state_cartesian = synodic2car(pre_DSM_state', pre_DSM_time, c.mu, theta_initial);
    post_DSM_state_cartesian = synodic2car(post_DSM_state', pre_DSM_time, c.mu, theta_initial);

    dv = norm(post_DSM_state_cartesian(4:6) - pre_DSM_state_cartesian(4:6))*(c.Vstar);

end
