function dv = objective_function_moon2comet_ga(x, c, target_pos_inertial, theta_inertial_final, q)
% objective_function_moon2comet_ga  GA cost for the Moon-to-comet leg:
% DSM delta-v magnitude plus a time-of-flight penalty (weight q).
%
% Inputs:
%   x                    - design variables [Vinf theta fpa oop tof beta]
%   c                    - constants struct
%   target_pos_inertial  - comet position, inertial [km]
%   theta_inertial_final - synodic frame angle at the comet epoch [rad]
%   q                    - weight on time of flight in the cost
%
% Outputs:
%   dv - scalar cost (DSM delta-v [km/s] + q*tof)

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

    theta_initial = theta_inertial_final - tof;
    pre_DSM_state_cartesian = synodic2car(pre_DSM_state', pre_DSM_time, c.mu, theta_initial);

    pre_DSM_state_cartesian_dimensionalized = [pre_DSM_state_cartesian(1:3) .* c.Lstar, pre_DSM_state_cartesian(4:6) .* c.Vstar];
    pre_DSM = pre_DSM_state_cartesian_dimensionalized;

    [V1, ~, ~, ~] = lambert(pre_DSM(1:3), target_pos_inertial, tof*(1-beta)*c.Tstar/(24*3600), 0, c.G*c.mSun);

    % dv = norm(V1 - pre_DSM(4:6)) + tof;
    dv = norm(V1 - pre_DSM(4:6)) + q*tof;

end
