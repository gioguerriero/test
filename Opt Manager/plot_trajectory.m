function plot_trajectory(out, S_halo, c, comet_name)
% plot_trajectory  Thin wrapper that plots the final full-ephemeris trajectory.
%
% Inputs:
%   out        - refined trajectory struct from run_refinement
%   S_halo     - halo orbit states
%   c          - constants struct
%   comet_name - string used for the figure title
%
% Outputs:
%   (none) - produces a figure

plot_full_trajectory(out, S_halo, c, comet_name);
end
