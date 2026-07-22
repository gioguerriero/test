function [initial_state, idx] = state_finder(theta_start, S_halo)
%======================================================================
% STATE_FINDER
%
% Selects the closest state on a Halo orbit based on a given angular
% position theta_start defined in the YZ plane.
%
% INPUTS:
%   theta_start  - Desired angular position in the YZ plane [rad]
%   S_halo       - Nx6 matrix containing the propagated Halo states
%                  over one full period (typically N ≈ 10000)
%
% OUTPUT:
%   initial_state - 1x6 state vector corresponding to the closest
%                   angular position in the YZ plane
%   idx           - index of the selected state within S_halo
%
% DESCRIPTION:
% - The angular coordinate is defined as:
%       theta = atan2(z, y)
% - Instead of computing atan2 for all states (computationally heavier),
%   we project each state onto the YZ plane and maximize alignment with
%   the desired unit direction.
% - This approach is efficient and well suited for optimization loops.
%
% ASSUMPTIONS:
% - The Halo orbit is ordered sequentially over one full period.
% - No state lies exactly at y = z = 0 (handled safely if it happens).
%======================================================================

    % Extract Y and Z coordinates
    y = S_halo(:,2);
    z = S_halo(:,3);

    % Desired unit direction in YZ plane
    uy_ref = cos(theta_start);
    uz_ref = sin(theta_start);

    % Compute radial distance in YZ plane
    r_yz = sqrt(y.^2 + z.^2);

    % Avoid division by zero
    r_yz(r_yz == 0) = eps;

    % Normalize YZ components
    uy = y ./ r_yz;
    uz = z ./ r_yz;

    % Dot product with reference direction
    % Equivalent to cos(angle_difference)
    dot_val = uy * uy_ref + uz * uz_ref;

    % Select index with maximum alignment
    [~, idx] = max(dot_val);

    % Return corresponding state
    initial_state = S_halo(idx,:)';

end


