function state = moon_state(theta, c)
% moon_state  Approximate Moon state in the Sun-Earth synodic frame as a
% function of its phase angle (circular Earth-Moon orbit model).
%
% Inputs:
%   theta - Moon phase angle in the synodic frame [rad]
%   c     - constants struct (uses rMoon_ad, mu)
%
% Outputs:
%   state - 1x6 Moon synodic state [x y z vx vy vz]

    % Moon position in the synodic frame
    x = c.rMoon_ad * cos(theta) + (1 - c.mu);
    y = c.rMoon_ad * sin(theta);
    z = 0;

    omega_moon = 365/27.3;   % Moon angular rate in synodic units (approx 13.37)
    vx = -omega_moon * c.rMoon_ad * sin(theta);
    vy =  omega_moon * c.rMoon_ad * cos(theta);
    vz = 0;

    state = [x; y; z; vx; vy; vz]';

end
