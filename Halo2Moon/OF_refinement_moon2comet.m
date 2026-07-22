function tof = OF_refinement_moon2comet(x, c)
% OF_refinement_moon2comet  Objective for the Halo-to-Moon refinement:
% total time of flight (Moon leg + Halo leg).
%
% Inputs:
%   x - design variables (uses x(4)=tof_moon, x(5)=tof_halo)
%   c - constants struct (unused)
%
% Outputs:
%   tof - total time of flight to minimize

tof_moon = x(4);
tof_halo = x(5);

tof = tof_moon + tof_halo;

end