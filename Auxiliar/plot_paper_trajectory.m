function plot_paper_trajectory(global_results, k, S_halo, c, selected_comet, unstable_dir, eps_vel_ms)
%PLOT_PAPER_TRAJECTORY  Publication-quality dual-panel CR3BP trajectory.
%
%  Reconstructs the CR3BP trajectory directly from global_results{k}.
%  Both panels are in the Sun-Earth synodic adimensional frame.
%
%  Layout:
%    Top panel    – full trajectory overview
%    Bottom panel – zoom on halo-orbit departure
%
%  Color convention (matches paper):
%    Halo orbit         – black, solid
%    Halo → Moon arcs  – blue
%    Moon → Comet arcs – green
%    DV vectors         – dark red quiver arrows
%    Moon orbit         – gray dashed thin
%
%  INPUTS
%    global_results – cell array from cr3bp_search
%    k              – index of the solution to plot
%    S_halo         – N×6 halo orbit states in synodic adim
%    c              – constants struct  (c.mu, c.Tstar, c.Lstar, c.Vstar, ...)
%    selected_comet – comet struct (.name, .comet_pos [km J2000], .epoch [ET s])
%    unstable_dir   – N×6 unstable manifold eigenvectors at each halo point
%    eps_vel_ms     – manifold injection magnitude [m/s]
%
%  USAGE (from main):
%    plot_paper_trajectory(global_results, 1, S_halo, c, selected_comet, ...
%                          unstable_dir, eps_vel_ms)
%
%  TUNING
%    DV_SCALE_OV  – arrow length in adim for overview panel
%    DV_SCALE_ZM  – arrow length in adim for zoom panel

%% ========================================================================
%  CONSTANTS AND STYLE
% =========================================================================
mu   = c.mu;
AU   = c.Lstar;   % 1 AU in km

COL_BLUE    = [0.18  0.42  0.75];
COL_GREEN   = [0.10  0.55  0.25];
COL_MOON_OR = [0.55  0.55  0.55];
COL_INJ     = [0.82  0.48  0.08];   % warm amber  – DV injection
COL_DSM1    = [0.55  0.20  0.62];   % muted purple – DSM1
COL_DSM2    = [0.78  0.10  0.12];   % dark red      – DSM2

LW_TRAJ  = 1.8;
LW_HALO  = 1.3;
LW_ORBIT = 0.7;

DV_SCALE_OV     = 0.007;   % arrow length in adim (overview, DV_inj and DSM1)
DV_SCALE_DSM2   = 0.06;    % arrow length in adim (overview, DSM2 – larger)
DV_SCALE_ZM     = 0.002;   % arrow length in adim (zoom panel)

eps_vel_ad = eps_vel_ms / (1e3 * c.Vstar);

%% ========================================================================
%  EXTRACT DESIGN VARIABLES FROM SOLUTION k
% =========================================================================
sol = global_results{k};

ref_fpa      = sol.x_opt(1);
ref_oop      = sol.x_opt(2);
ref_halo_th  = sol.x_opt(3);
ref_tof_moon = sol.x_opt(4);
ref_tof_halo = sol.x_opt(5);

rf         = sol.result_fmincon;
Vinf_m2c   = rf(1);
theta_m2c  = rf(2);
fpa_m2c    = rf(3);
oop_m2c    = rf(4);
tof_m2c_ad = rf(5);
beta_m2c   = rf(6);
dv_m2c_ad  = rf(7:9)';

fvr = sol.fmincon_vinf_results;

%% ========================================================================
%  CR3BP ARC PROPAGATION (synodic adim)
% =========================================================================
opt_plot = odeset('AbsTol', 1e-10, 'RelTol', 1e-10);
N_arc    = 500;

% Arc 1: Halo departure → DSM1
[init_st, idx_h] = state_finder(ref_halo_th, S_halo);
vu  = unstable_dir(idx_h, :)';
vu  = vu / norm(vu);
S0_arc1  = [init_st(1:3); init_st(4:6) + eps_vel_ad * vu];
[~, S1]  = ode45(@(t,S) CR3BP(t, S, mu), linspace(0, ref_tof_halo, N_arc), S0_arc1, opt_plot);

% Arc 2: DSM1 → Moon flyby (backwards from Moon, then reversed)
moon_syn_h2m  = moon_state(fvr.theta_moon, c)';
v_inf_h2m     = vinf_rotation(moon_syn_h2m, fvr.vinf_fmincon, ref_fpa, ref_oop);
S0_flyby      = moon_syn_h2m;
S0_flyby(4:6) = S0_flyby(4:6) + v_inf_h2m ./ c.Vstar;
[~, S2b]      = ode45(@(t,S) CR3BP(t, S, mu), linspace(0, -ref_tof_moon, N_arc), S0_flyby, opt_plot);
S2            = flipud(S2b);

% Arc 3: Moon flyby → DSM2
moon_syn_m2c = moon_state(theta_m2c, c);
vinf_m2c_vec = vinf_rotation(moon_syn_m2c', Vinf_m2c, fpa_m2c, oop_m2c);
S0_arc3      = [moon_syn_m2c(1:3)'; moon_syn_m2c(4:6)' + vinf_m2c_vec ./ c.Vstar];
[~, S3]      = ode45(@(t,S) CR3BP(t, S, mu), linspace(0, tof_m2c_ad*beta_m2c, N_arc), S0_arc3, opt_plot);

% Arc 4: DSM2 → Comet
S0_arc4 = [S3(end,1:3)'; S3(end,4:6)' + dv_m2c_ad];
[~, S4] = ode45(@(t,S) CR3BP(t, S, mu), linspace(0, tof_m2c_ad*(1-beta_m2c), N_arc), S0_arc4, opt_plot);

%% ========================================================================
%  MANEUVER VECTORS (synodic adim)
% =========================================================================
dv_inj_syn  = eps_vel_ad * vu;
dv_dsm1_syn = S2(1,4:6)'  - S1(end,4:6)';
dv_dsm2_syn = dv_m2c_ad;

dv_inj_ms  = eps_vel_ms;
dv_dsm1_ms = norm(dv_dsm1_syn) * c.Vstar * 1e3;
dv_dsm2_ms = norm(dv_dsm2_syn) * c.Vstar * 1e3;

% Maneuver positions (synodic)
pos_dep_syn  = S0_arc1(1:3)';
pos_dsm1_syn = S1(end, 1:3);
pos_dsm2_syn = S3(end, 1:3);

% Unit directions for arrows
dir_inj_syn  = (dv_inj_syn  / norm(dv_inj_syn))';
dir_dsm1_syn = (dv_dsm1_syn / norm(dv_dsm1_syn))';
dir_dsm2_syn = (dv_dsm2_syn / norm(dv_dsm2_syn))';

%% ========================================================================
%  REFERENCE GEOMETRY (synodic adim)
% =========================================================================
% Lagrange points
rL_ad  = (mu/3)^(1/3);
L1_syn = [1-mu-rL_ad, 0, 0];
L2_syn = [1-mu+rL_ad, 0, 0];

% Moon orbit: approximate circle around Earth in synodic
r_moon_adim   = 384400 / AU;
theta_circle  = linspace(0, 2*pi, 200);
moon_orb_x    = (1-mu) + r_moon_adim * cos(theta_circle);
moon_orb_y    =           r_moon_adim * sin(theta_circle);

% Halo orbit (N×6 synodic)
if size(S_halo, 2) == 6,  S_halo_syn = S_halo;
else,                      S_halo_syn = S_halo';
end

% Comet position in synodic adim (rotate J2000 position to synodic at arrival epoch)
earth_arr  = cspice_spkezr('EARTH', selected_comet.epoch, 'ECLIPJ2000', 'NONE', 'SUN');
theta_arr  = atan2(earth_arr(2), earth_arr(1));
Lstar_arr  = norm(earth_arr(1:2));
Rz         = [ cos(theta_arr)  sin(theta_arr)  0;
              -sin(theta_arr)  cos(theta_arr)  0;
               0               0               1];
comet_syn  = Rz * selected_comet.comet_pos(:) / Lstar_arr;  % 3×1

% Departure epoch (approximate, from arrival backward)
tof_total_s = (ref_tof_halo + ref_tof_moon + tof_m2c_ad) * c.Tstar;
epoch_dep   = selected_comet.epoch - tof_total_s;

%% ========================================================================
%  FIGURE LAYOUT  (overview top, zoom bottom)
% =========================================================================
hfig = figure('Name', ...
    sprintf('Trajectory %s – solution %d', selected_comet.name, k), ...
    'Color', 'w', 'Position', [100, 50, 900, 1050]);

ax1 = axes('Parent', hfig, 'Position', [0.08, 0.50, 0.88, 0.46]);  % top (slightly larger)
ax2 = axes('Parent', hfig, 'Position', [0.08, 0.05, 0.88, 0.42]);  % bottom

%% ========================================================================
%  PANEL 1 (top) – FULL TRAJECTORY OVERVIEW  (synodic adim)
% =========================================================================
hold(ax1, 'on');  grid(ax1, 'on');  axis(ax1, 'equal');
ax1.GridAlpha = 0.20;

% Halo orbit – black (not in legend at this scale)
plot3(ax1, S_halo_syn(:,1), S_halo_syn(:,2), S_halo_syn(:,3), ...
    'k-', 'LineWidth', LW_HALO, 'HandleVisibility', 'off');

% Moon orbit – gray dashed (not in legend at this scale)
plot3(ax1, moon_orb_x, moon_orb_y, zeros(size(theta_circle)), ...
    '--', 'Color', COL_MOON_OR, 'LineWidth', LW_ORBIT, 'HandleVisibility', 'off');

% Halo → Moon – blue (arcs 1 and 2)
plot3(ax1, S1(:,1), S1(:,2), S1(:,3), ...
    '-', 'Color', COL_BLUE, 'LineWidth', LW_TRAJ, 'DisplayName', 'Halo \rightarrow Moon');
plot3(ax1, S2(:,1), S2(:,2), S2(:,3), ...
    '-', 'Color', COL_BLUE, 'LineWidth', LW_TRAJ, 'HandleVisibility', 'off');

% Moon → Comet – green (arcs 3 and 4)
plot3(ax1, S3(:,1), S3(:,2), S3(:,3), ...
    '-', 'Color', COL_GREEN, 'LineWidth', LW_TRAJ, 'DisplayName', 'Moon \rightarrow Comet');
plot3(ax1, S4(:,1), S4(:,2), S4(:,3), ...
    '-', 'Color', COL_GREEN, 'LineWidth', LW_TRAJ, 'HandleVisibility', 'off');

% DV injection arrow (not in legend at this scale)
quiver3(ax1, pos_dep_syn(1),  pos_dep_syn(2),  pos_dep_syn(3), ...
    dir_inj_syn(1)*DV_SCALE_OV,  dir_inj_syn(2)*DV_SCALE_OV,  dir_inj_syn(3)*DV_SCALE_OV, ...
    0, 'Color', COL_INJ, 'LineWidth', 2.0, 'MaxHeadSize', 0.6, 'HandleVisibility', 'off');

% DSM1 arrow (not in legend at this scale)
quiver3(ax1, pos_dsm1_syn(1), pos_dsm1_syn(2), pos_dsm1_syn(3), ...
    dir_dsm1_syn(1)*DV_SCALE_OV, dir_dsm1_syn(2)*DV_SCALE_OV, dir_dsm1_syn(3)*DV_SCALE_OV, ...
    0, 'Color', COL_DSM1, 'LineWidth', 2.0, 'MaxHeadSize', 0.6, 'HandleVisibility', 'off');

% DSM2 arrow (larger than DV_inj and DSM1)
quiver3(ax1, pos_dsm2_syn(1), pos_dsm2_syn(2), pos_dsm2_syn(3), ...
    dir_dsm2_syn(1)*DV_SCALE_DSM2, dir_dsm2_syn(2)*DV_SCALE_DSM2, dir_dsm2_syn(3)*DV_SCALE_DSM2, ...
    0, 'Color', COL_DSM2, 'LineWidth', 2.5, 'MaxHeadSize', 0.4, 'DisplayName', 'DSM2');

% Sun
plot3(ax1, -mu, 0, 0, 'o', 'Color', [1.00 0.75 0.00], ...
    'MarkerFaceColor', [1.00 0.75 0.00], 'MarkerSize', 14, 'DisplayName', 'Sun');

% Earth
plot3(ax1, 1-mu, 0, 0, 'o', 'Color', [0.20 0.50 1.00], ...
    'MarkerFaceColor', [0.20 0.50 1.00], 'MarkerSize', 8, 'DisplayName', 'Earth');

% L1, L2 (not in legend at this scale)
plot3(ax1, L1_syn(1), L1_syn(2), L1_syn(3), 'k+', 'MarkerSize', 9, ...
    'LineWidth', 1.5, 'HandleVisibility', 'off');
plot3(ax1, L2_syn(1), L2_syn(2), L2_syn(3), 'kx', 'MarkerSize', 9, ...
    'LineWidth', 1.5, 'HandleVisibility', 'off');

% Comet
plot3(ax1, comet_syn(1), comet_syn(2), comet_syn(3), 'p', ...
    'Color', [0.55 0.10 0.60], 'MarkerFaceColor', [0.55 0.10 0.60], ...
    'MarkerSize', 13, 'DisplayName', selected_comet.name);

xlabel(ax1, 'x  [adim]');  ylabel(ax1, 'y  [adim]');
legend(ax1, 'show', 'Location', 'best', 'FontSize', 9, 'Interpreter', 'tex');
view(ax1, 0, 90);

%% ========================================================================
%  PANEL 2 (bottom) – ZOOM ON HALO DEPARTURE  (synodic adim)
% =========================================================================
hold(ax2, 'on');  grid(ax2, 'on');  axis(ax2, 'equal');
ax2.GridAlpha = 0.20;

% Same content as overview — xlim/ylim crop to the halo region

% Halo orbit – black
plot3(ax2, S_halo_syn(:,1), S_halo_syn(:,2), S_halo_syn(:,3), ...
    'k-', 'LineWidth', LW_HALO, 'DisplayName', 'Halo orbit');

% Moon orbit – gray dashed
plot3(ax2, moon_orb_x, moon_orb_y, zeros(size(theta_circle)), ...
    '--', 'Color', COL_MOON_OR, 'LineWidth', LW_ORBIT, 'DisplayName', 'Moon orbit');

% Halo → Moon – blue (arcs 1 and 2)
plot3(ax2, S1(:,1), S1(:,2), S1(:,3), ...
    '-', 'Color', COL_BLUE, 'LineWidth', LW_TRAJ, 'DisplayName', 'Halo \rightarrow Moon');
plot3(ax2, S2(:,1), S2(:,2), S2(:,3), ...
    '-', 'Color', COL_BLUE, 'LineWidth', LW_TRAJ, 'HandleVisibility', 'off');

% Moon → Comet – green (arcs 3 and 4)
plot3(ax2, S3(:,1), S3(:,2), S3(:,3), ...
    '-', 'Color', COL_GREEN, 'LineWidth', LW_TRAJ, 'DisplayName', 'Moon \rightarrow Comet');
plot3(ax2, S4(:,1), S4(:,2), S4(:,3), ...
    '-', 'Color', COL_GREEN, 'LineWidth', LW_TRAJ, 'HandleVisibility', 'off');

% DV injection arrow
quiver3(ax2, pos_dep_syn(1), pos_dep_syn(2), pos_dep_syn(3), ...
    dir_inj_syn(1)*DV_SCALE_ZM, dir_inj_syn(2)*DV_SCALE_ZM, dir_inj_syn(3)*DV_SCALE_ZM, ...
    0, 'Color', COL_INJ, 'LineWidth', 2.5, 'MaxHeadSize', 0.8, ...
    'DisplayName', '\DeltaV_{inj}');

% DSM1 arrow
quiver3(ax2, pos_dsm1_syn(1), pos_dsm1_syn(2), pos_dsm1_syn(3), ...
    dir_dsm1_syn(1)*DV_SCALE_ZM, dir_dsm1_syn(2)*DV_SCALE_ZM, dir_dsm1_syn(3)*DV_SCALE_ZM, ...
    0, 'Color', COL_DSM1, 'LineWidth', 2.5, 'MaxHeadSize', 0.8, ...
    'DisplayName', 'DSM1');

% Departure point marker
plot3(ax2, pos_dep_syn(1), pos_dep_syn(2), pos_dep_syn(3), 's', ...
    'Color', COL_BLUE, 'MarkerFaceColor', COL_BLUE, 'MarkerSize', 7, ...
    'HandleVisibility', 'off');

% Moon flyby marker
plot3(ax2, S2(end,1), S2(end,2), S2(end,3), 'o', ...
    'Color', COL_MOON_OR, 'MarkerFaceColor', COL_MOON_OR, ...
    'MarkerSize', 9, 'DisplayName', 'Moon flyby');

% Sun
plot3(ax2, -mu, 0, 0, 'o', 'Color', [1.00 0.75 0.00], ...
    'MarkerFaceColor', [1.00 0.75 0.00], 'MarkerSize', 14, 'HandleVisibility', 'off');

% Earth
plot3(ax2, 1-mu, 0, 0, 'o', 'Color', [0.20 0.50 1.00], ...
    'MarkerFaceColor', [0.20 0.50 1.00], 'MarkerSize', 8, 'DisplayName', 'Earth');

% L1, L2
plot3(ax2, L1_syn(1), L1_syn(2), L1_syn(3), 'k+', 'MarkerSize', 9, ...
    'LineWidth', 1.5, 'DisplayName', 'L_1');
plot3(ax2, L2_syn(1), L2_syn(2), L2_syn(3), 'kx', 'MarkerSize', 9, ...
    'LineWidth', 1.5, 'DisplayName', 'L_2');

% Zoom window: wider in x to show arc departing from halo, tight in y
margin_x = 0.04;   % horizontal extent to follow the departing arc
margin_y = 0.015;  % tight vertical margin
xlim(ax2, [min(S_halo_syn(:,1)) - margin_x,  max(S_halo_syn(:,1)) + margin_x]);
ylim(ax2, [min(S_halo_syn(:,2)) - margin_y,  max(S_halo_syn(:,2)) + margin_y]);

xlabel(ax2, 'x  [adim]');  ylabel(ax2, 'y  [adim]');
legend(ax2, 'show', 'Location', 'best', 'FontSize', 9, 'Interpreter', 'tex');
view(ax2, 0, 90);

end
