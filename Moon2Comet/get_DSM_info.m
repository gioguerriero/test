function dv = get_DSM_info(result, c, theta_inertial_final, target_pos_inertial)
% get_DSM_info  Compute the deep-space maneuver (DSM) delta-v for a Moon-to-comet
% leg by propagating to the DSM point and solving Lambert to the comet.
%
% Inputs:
%   result               - optimized M2C variables [Vinf theta fpa oop tof beta]
%   c                    - constants struct
%   theta_inertial_final - synodic frame angle at the comet epoch [rad]
%   target_pos_inertial  - comet position, inertial [km]
%
% Outputs:
%   dv - DSM delta-v in the synodic (non-dimensional) frame [1x3]

    Vinf  = result(1);
    theta = result(2);
    fpa = result(3);
    out_of_plane  = result(4);
    tof   = result(5);
    beta  = result(6);

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

    % Post-maneuver state (inertial)
    post_DSM_state_intertial = [pre_DSM(1:3)./c.Lstar, V1./c.Vstar];

    % Convert the post-maneuver state to the synodic frame
    post_DSM_state_synodic = car2synodic(post_DSM_state_intertial, pre_DSM_time, c.mu, theta_initial);

    % Delta-v in the synodic (non-dimensional) frame
    dv = post_DSM_state_synodic(4:6) - pre_DSM_state(4:6)';

end