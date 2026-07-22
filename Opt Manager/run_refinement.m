function out = run_refinement(c, best, selected_comet, S_halo, unstable_dir, rp)
%RUN_REFINEMENT  Full ephemeris refinement of a CR3BP solution.
%
%  Inputs:
%    c              - constants struct
%    best           - solution struct from cr3bp_search / pareto_data
%    selected_comet - comet struct
%    S_halo         - halo orbit states
%    unstable_dir   - unstable eigenvectors
%    rp             - refine_params struct (see main.m)
%
%  Output:
%    out  - fully refined trajectory struct from variables_organizer_refined

fprintf('\n==========================================================\n');
fprintf('  REFINEMENT: %s | TOF=%.2f d | dv_max=%.0f m/s | method=%s\n', ...
    selected_comet.name, best.tof_total_days, best.maximum_dv*1e3, ...
    choose_str(rp.use_fast, 'fast', 'global'));
fprintf('==========================================================\n');

eps_vel_ms = rp.eps_vel_ms;
eps_vel_ad = eps_vel_ms / (1e3 * c.Vstar);

%% ================================================================
%  ARC PROPAGATION (CR3BP, synodic adim)
% ================================================================
opt_traj = odeset('AbsTol', 1e-10, 'RelTol', 1e-10);
N_arc    = 2000;
N_full   = 5000;

% H2M design variables
ref_fpa      = best.x_opt(1);
ref_oop      = best.x_opt(2);
ref_halo_th  = best.x_opt(3);
ref_tof_moon = best.x_opt(4);
ref_tof_halo = best.x_opt(5);

% M2C design variables
rf         = best.result_fmincon;
Vinf_m2c   = rf(1);
theta_m2c  = rf(2);
fpa_m2c    = rf(3);
oop_m2c    = rf(4);
tof_m2c_ad = rf(5);
beta_m2c   = rf(6);
dv_m2c_ad  = rf(7:9)';

fvr = best.fmincon_vinf_results;

% Arc 1: Halo departure -> DSM1
[init_st, idx_h] = state_finder(ref_halo_th, S_halo);
vu = unstable_dir(idx_h,:)';
vu = vu / norm(vu);

S_pre_inj_syn = [init_st(1:3); init_st(4:6)];
S0_arc1       = [init_st(1:3); init_st(4:6) + eps_vel_ad * vu];
[t1, S1]      = ode45(@(t,S) CR3BP(t,S,c.mu), linspace(0, ref_tof_halo, N_arc), S0_arc1, opt_traj);

% Arc 2: DSM1 -> Moon flyby (back-prop then reversed)
moon_syn  = moon_state(fvr.theta_moon, c)';
v0_h2m    = vinf_rotation(moon_syn, fvr.vinf_fmincon, fvr.fpa_fmincon, fvr.out_of_plane_fmincon);
v_inf_h2m = vinf_rotation(moon_syn, norm(v0_h2m), ref_fpa, ref_oop);
S0_flyby  = moon_syn;
S0_flyby(4:6) = S0_flyby(4:6) + v_inf_h2m ./ c.Vstar;

[~, S2b] = ode45(@(t,S) CR3BP(t,S,c.mu), linspace(0, -ref_tof_moon, N_arc), S0_flyby, opt_traj);
S2        = flipud(S2b);
t2_local  = linspace(0, ref_tof_moon, N_arc)';
t2        = t2_local + ref_tof_halo;

% Arc 3: Moon flyby -> DSM2
moon_syn_m2c = moon_state(theta_m2c, c);
vinf_vec_m2c = vinf_rotation(moon_syn_m2c', Vinf_m2c, fpa_m2c, oop_m2c);
S0_arc3      = [moon_syn_m2c(1:3)'; moon_syn_m2c(4:6)' + vinf_vec_m2c ./ c.Vstar];
[t3_local, S3] = ode45(@(t,S) CR3BP(t,S,c.mu), linspace(0, tof_m2c_ad*beta_m2c, N_arc), S0_arc3, opt_traj);
t3 = t3_local + ref_tof_halo + ref_tof_moon;

% Arc 4: DSM2 -> Comet
S0_arc4 = [S3(end,1:3)'; S3(end,4:6)' + dv_m2c_ad];
[t4_local, S4] = ode45(@(t,S) CR3BP(t,S,c.mu), linspace(0, tof_m2c_ad*(1-beta_m2c), N_arc), S0_arc4, opt_traj);
t4 = t4_local + ref_tof_halo + ref_tof_moon + tof_m2c_ad*beta_m2c;

% Maneuver vectors (synodic adim)
dv_inj_ad_syn   = eps_vel_ad * vu;
dv_dsm1_ad_syn  = S2(1,4:6)' - S1(end,4:6)';
dv_flyby_ad_syn = S0_arc3(4:6) - S2(end,4:6)';
dv_dsm2_ad_syn  = dv_m2c_ad;

S_pre_inj_syn_6  = S_pre_inj_syn;
S_pre_dsm1_syn   = S1(end,:)';
S_pre_flyby_syn  = S2(end,:)';
S_pre_dsm2_syn   = S3(end,:)';

% Concatenate arcs and resample
t_cat       = [t1; t2(2:end); t3(2:end); t4(2:end)];
S_cat       = [S1; S2(2:end,:); S3(2:end,:); S4(2:end,:)];
t_cat_s     = t_cat * c.Tstar;
tof_total_s = t_cat_s(end);

t_traj_s   = linspace(0, tof_total_s, N_full)';
S_traj_syn = zeros(N_full, 6);
for kk = 1:6
    S_traj_syn(:,kk) = interp1(t_cat_s, S_cat(:,kk), t_traj_s, 'pchip');
end

%% ================================================================
%  EPOCHS (ET seconds)
% ================================================================
epoch_arr   = selected_comet.epoch;
epoch_vec   = epoch_arr - (tof_total_s - t_traj_s);
epoch_dep   = epoch_vec(1);
epoch_dsm1  = epoch_dep   + ref_tof_halo          * c.Tstar;
epoch_flyby = epoch_dsm1  + ref_tof_moon           * c.Tstar;
epoch_dsm2  = epoch_flyby + tof_m2c_ad * beta_m2c  * c.Tstar;

%% ================================================================
%  INTERMEDIATE NODES (for multiple shooting initial guess)
% ================================================================
n_h_dsm1 = 3; n_dsm1_f = 3; n_f_dsm2 = 5; n_dsm2_c = 5;

t_nodes_h_dsm1 = (1:n_h_dsm1) * (ref_tof_halo            / (n_h_dsm1+1));
t_nodes_dsm1_f = (1:n_dsm1_f) * (ref_tof_moon             / (n_dsm1_f+1));
t_nodes_f_dsm2 = (1:n_f_dsm2) * (tof_m2c_ad*beta_m2c      / (n_f_dsm2+1));
t_nodes_dsm2_c = (1:n_dsm2_c) * (tof_m2c_ad*(1-beta_m2c)  / (n_dsm2_c+1));

S_nodes_h_dsm1_syn = zeros(n_h_dsm1, 6);
S_nodes_dsm1_f_syn = zeros(n_dsm1_f, 6);
S_nodes_f_dsm2_syn = zeros(n_f_dsm2, 6);
S_nodes_dsm2_c_syn = zeros(n_dsm2_c, 6);
for kk = 1:6
    S_nodes_h_dsm1_syn(:,kk) = interp1(t1,       S1(:,kk), t_nodes_h_dsm1, 'pchip');
    S_nodes_dsm1_f_syn(:,kk) = interp1(t2_local, S2(:,kk), t_nodes_dsm1_f, 'pchip');
    S_nodes_f_dsm2_syn(:,kk) = interp1(t3_local, S3(:,kk), t_nodes_f_dsm2, 'pchip');
    S_nodes_dsm2_c_syn(:,kk) = interp1(t4_local, S4(:,kk), t_nodes_dsm2_c, 'pchip');
end

epoch_nodes_h_dsm1 = epoch_dep   + t_nodes_h_dsm1 * c.Tstar;
epoch_nodes_dsm1_f = epoch_dsm1  + t_nodes_dsm1_f * c.Tstar;
epoch_nodes_f_dsm2 = epoch_flyby + t_nodes_f_dsm2 * c.Tstar;
epoch_nodes_dsm2_c = epoch_dsm2  + t_nodes_dsm2_c * c.Tstar;

%% ================================================================
%  SYNODIC -> SUN-CENTERED J2000 (batch SPICE call)
% ================================================================
epoch_all = [epoch_vec(:);
             epoch_dep; epoch_dsm1; epoch_flyby; epoch_dsm2;
             epoch_nodes_h_dsm1(:);
             epoch_nodes_dsm1_f(:);
             epoch_nodes_f_dsm2(:);
             epoch_nodes_dsm2_c(:)];

S_all_syn = [S_traj_syn;
             S_pre_inj_syn_6.';
             S_pre_dsm1_syn.';
             S_pre_flyby_syn.';
             S_pre_dsm2_syn.';
             S_nodes_h_dsm1_syn;
             S_nodes_dsm1_f_syn;
             S_nodes_f_dsm2_syn;
             S_nodes_dsm2_c_syn];

[S_all_J2000, scales] = synodic2sun_J2000(S_all_syn, epoch_all, c.mu);

idx_traj     = 1:N_full;
idx_ev_inj   = N_full+1;  idx_ev_dsm1  = N_full+2;
idx_ev_flyby = N_full+3;  idx_ev_dsm2  = N_full+4;
nh_s = N_full+5;          nh_e = nh_s+n_h_dsm1-1;
df_s = nh_e+1;            df_e = df_s+n_dsm1_f-1;
fd_s = df_e+1;            fd_e = fd_s+n_f_dsm2-1;
dc_s = fd_e+1;            dc_e = dc_s+n_dsm2_c-1;

S_traj            = S_all_J2000(idx_traj,   :);
S_pre_inj_J2000   = S_all_J2000(idx_ev_inj,  :).';
S_pre_dsm1_J2000  = S_all_J2000(idx_ev_dsm1, :).';
S_pre_flyby_J2000 = S_all_J2000(idx_ev_flyby,:).';
S_pre_dsm2_J2000  = S_all_J2000(idx_ev_dsm2, :).';

% DV conversion: synodic adim -> J2000 km/s
dv_inj_vec_kms   = scales.R(:,:,idx_ev_inj)   * (dv_inj_ad_syn   * scales.V(idx_ev_inj));
dv_dsm1_vec_kms  = scales.R(:,:,idx_ev_dsm1)  * (dv_dsm1_ad_syn  * scales.V(idx_ev_dsm1));
dv_flyby_vec_kms = scales.R(:,:,idx_ev_flyby) * (dv_flyby_ad_syn * scales.V(idx_ev_flyby));
dv_dsm2_vec_kms  = scales.R(:,:,idx_ev_dsm2)  * (dv_dsm2_ad_syn  * scales.V(idx_ev_dsm2));

% TOF segments [days]
tof_halo2dsm1_d  = ref_tof_halo           * c.Tstar / 86400;
tof_dsm12flyby_d = ref_tof_moon           * c.Tstar / 86400;
tof_flyby2dsm2_d = tof_m2c_ad*beta_m2c   * c.Tstar / 86400;
tof_dsm22comet_d = tof_m2c_ad*(1-beta_m2c)* c.Tstar / 86400;

%% Print CR3BP summary
fprintf('\n--- TOF segments ---\n');
fprintf('  Total:           %.4f d\n', tof_total_s/86400);
fprintf('  Halo  -> DSM1:   %.4f d\n', tof_halo2dsm1_d);
fprintf('  DSM1  -> Flyby:  %.4f d\n', tof_dsm12flyby_d);
fprintf('  Flyby -> DSM2:   %.4f d\n', tof_flyby2dsm2_d);
fprintf('  DSM2  -> Comet:  %.4f d\n', tof_dsm22comet_d);
fprintf('\n--- DV [km/s] ---\n');
fprintf('  Injection: %.5f\n', norm(dv_inj_vec_kms));
fprintf('  DSM1:      %.5f\n', norm(dv_dsm1_vec_kms));
fprintf('  Flyby:     %.5f  (fictitious)\n', norm(dv_flyby_vec_kms));
fprintf('  DSM2:      %.5f\n', norm(dv_dsm2_vec_kms));

% Flyby bending angle and periapsis altitude
vinf_mag  = norm(v0_h2m);
cos_delta = max(-1, min(1, dot(v0_h2m, v_inf_h2m) / vinf_mag^2));
delta_rad = acos(cos_delta);
rp_flyby  = c.muMoon / vinf_mag^2 * (1/sin(delta_rad/2) - 1);
fprintf('\n--- Lunar flyby (zeroSOI) ---\n');
fprintf('  |v_inf| = %.4f km/s | delta = %.4f deg | h_p = %.2f km\n', ...
    vinf_mag, rad2deg(delta_rad), rp_flyby - c.rMoon);

% Heliocentric plot (CR3BP-based)
plot_trajectory_heliocentric(S_traj, t_traj_s, epoch_dep, epoch_dsm1, ...
    epoch_flyby, epoch_dsm2, S_pre_inj_J2000, S_pre_dsm1_J2000, ...
    S_pre_flyby_J2000, S_pre_dsm2_J2000, selected_comet);

%% ================================================================
%  EPHEMERIS OPTIMIZATION
% ================================================================
S_pre_inj_manual  = S_pre_inj_syn_6;
target_position   = selected_comet.comet_pos;
epoch_comet_flyby = selected_comet.epoch;
multiple_shooting = rp.multiple_shooting;
k_vec             = rp.k_vec;

if multiple_shooting
    t_halo   = 0.5;   % [days]
    d_max    = 111400000 / c.Lstar;
    dt_max   = 3;
    % Epochs expressed as "days before comet arrival": x = (epoch_comet - epoch)/86400
    % Base (17): the dep/dsm1/flyby/dsm2 event epochs are now ALL explicit;
    % the node epochs are no longer design variables.
    x0_ephe  = [dv_inj_vec_kms.*1e3; dv_dsm1_vec_kms.*10; ...
                (epoch_comet_flyby - epoch_flyby)/86400; S_traj(end,4:6)'; dv_dsm2_vec_kms.*10; ...
                t_halo; (epoch_comet_flyby - epoch_dep)/86400; ...
                (epoch_comet_flyby - epoch_dsm1)/86400; (epoch_comet_flyby - epoch_dsm2)/86400];


    % --- flyby phase correction at CONSTANT total TOF ---
    % Advance the flyby by d_shift days: DSM1->Flyby shortens, Flyby->DSM2
    % lengthens, departure and comet stay fixed (total TOF unchanged).
    d_shift = 0;   % [days]  (was a hardcoded manual offset)
    x0_ephe(7)  = x0_ephe(7) + d_shift;        % only the flyby, as a design variable
    epoch_flyby = epoch_flyby - d_shift*86400; % consistency: same shift on the scalar [s]
    % NB: do NOT touch epoch_dep, epoch_dsm1, epoch_dsm2 -> total TOF constant

    [x0_complete, f_nodes] = build_initial_guess_MS(x0_ephe, k_vec, S_traj_syn, epoch_vec, ...
        S_pre_dsm1_syn, S_pre_dsm2_syn, epoch_dep, epoch_dsm1, ...
        epoch_flyby, epoch_dsm2, epoch_comet_flyby);

    % Plot to visualise how the trajectory is split and check the chosen k_vec
    plot_from_state_vector_MS(x0_complete, k_vec, f_nodes, S_traj_syn, c, target_position, epoch_comet_flyby, S_halo);

    optimality_tolerance = 1e-3;
    constraint_tolerance = 1e-3;

    [result_ephe, x_opt] = optimization_ephe(x0_complete, S_pre_inj_manual, rp.max_dv, ...
        target_position, epoch_comet_flyby, rp.hp, c, rp.max_dv_inj, ...
        multiple_shooting, k_vec, f_nodes, d_max, dt_max, optimality_tolerance, constraint_tolerance, rp.min_days_between);

    % Either reduce the nodes, or take the roughly-correct guess and zero the
    % cost function so it only closes the gaps, hoping the delta-v stays put

else
    t_halo   = 0.5;
    x0_ephe  = [dv_inj_vec_kms.*1e3; tof_halo2dsm1_d; dv_dsm1_vec_kms.*10; ...
                tof_dsm12flyby_d; epoch_flyby/1e8; S_traj(end,4:6)'; ...
                tof_dsm22comet_d/100; dv_dsm2_vec_kms.*10; t_halo];
    [result_ephe, ~] = optimization_ephe(x0_ephe, S_pre_inj_manual, rp.max_dv, ...
        target_position, epoch_comet_flyby, rp.hp, c, rp.max_dv_inj, ...
        multiple_shooting, [], [], []);
    x_opt = result_ephe.bestfeasible.x;
end

out_cr3bp = variables_organizer(x_opt, k_vec, f_nodes, S_pre_inj_manual, ...
   target_position, epoch_comet_flyby, c);
save("out_cr3bp.mat", "out_cr3bp")

info = verify_flyby_geometry(out_cr3bp, c);
write_python_inputs(out_cr3bp, c)
write_python_inputs_serot(out_cr3bp, c) 
[h_min_earth, h_flyby, info] = check_altitude(out_cr3bp, c, 'S_halo', S_halo);

plot_from_state_vector_MS(x_opt, k_vec, f_nodes, S_traj_syn, c, target_position, epoch_comet_flyby, S_halo);


%% ================================================================
%  B-PLANE TARGETING
% ================================================================
[~, BT_target, BR_target, ~, ~, ~, ~, ~] = bplane_from_vinf( ...
    out_cr3bp.flyby.vinf_in, out_cr3bp.flyby.vinf_out, c.muMoon);

t_tcm1              = rp.t_tcm1;
epoch_flyby_ref     = out_cr3bp.flyby.epoch;
epoch_tcm1          = epoch_flyby_ref - t_tcm1 * 86400;

opt_prop = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
[~, S_back] = ode45(@(t,s) NBODY_J2000(t, s, epoch_flyby_ref, c), ...
    [0, -t_tcm1*86400], out_cr3bp.flyby.state_pre.J2000(:), opt_prop);
state_pre_tcm1 = S_back(end,:)';

[dv_tcm, ~, ~] = bplane_tcm(state_pre_tcm1, epoch_tcm1, epoch_flyby_ref, ...
    'MOON', BT_target, BR_target, c.muMoon, c);

%% ================================================================
%  POST-FLYBY REFINEMENT
% ================================================================
if rp.use_fast
    % Fast path: refinement_post_flyby
    t_tcm2      = 5;   % [days] post flyby
    tcm2_guess  = [0.001; 0.001; 0.001];
    x0_post     = [tcm2_guess; out_cr3bp.dsm2.dv.J2000_kms; out_cr3bp.dsm2.epoch/1e8];

    [~, x_post_flyby, info_post] = refinement_post_flyby(x0_post, c, state_pre_tcm1, dv_tcm, ...
        epoch_tcm1, target_position, epoch_comet_flyby, t_tcm2, epoch_flyby_ref);

else
    % Slow path: refinement_global_post_flyby + MBH
    t_tcm2 = 3;   % [days] post flyby                               !!!!!!!!!!!!!!!! note: was 10 days !!!!!!
    [~, S_post] = ode45(@(t,s) NBODY_J2000(t, s, epoch_flyby_ref, c), ...
        [0, t_tcm2*86400], out_cr3bp.flyby.state_post.J2000(:), opt_prop);
    S_post_tcm2 = S_post(end,:)';

    % x0_global = [BT_target/1e3; BR_target/1e3; out_cr3bp.dsm2.dv.J2000_kms; ...
    %             out_cr3bp.dsm2.epoch/1e8; S_post_tcm2(1:3)./1e8; S_post_tcm2(4:6)./10];

    tcm2_guess = [0.001; 0.001; 0.001];
    dsm2_epoch_variation = 0.1; %[in days]
    x0_global = [BT_target/5000; BR_target/5000; out_cr3bp.dsm2.dv.J2000_kms./0.1; ...
                dsm2_epoch_variation; out_cr3bp.comet.state.J2000(4:6)'./10; tcm2_guess./0.1];

    [mbh_out, x_post_flyby, info_post] = refinement_global_post_flyby(x0_global, c, epoch_tcm1, ...
        target_position, epoch_comet_flyby, t_tcm2, epoch_flyby_ref, state_pre_tcm1, rp.hp, out_cr3bp.dsm2.epoch);

end

%% ================================================================
%  PRE-FLYBY REFINEMENT
% ================================================================
dsm1_epoch_variation = 0.1; %[in days]
x0_pre = [out_cr3bp.injection.J2000_kms; ...
          out_cr3bp.dsm1.dv.J2000_kms;   ...
          dsm1_epoch_variation];

[~, x_pre_flyby, info_pre] = refinement_pre_flyby(x0_pre, c, state_pre_tcm1, epoch_tcm1, ...
    out_cr3bp.departure.state.J2000, out_cr3bp.departure.epoch, out_cr3bp.dsm1.epoch);

%% ================================================================
%  FINAL ASSEMBLY
% ================================================================
out = variables_organizer_refined(out_cr3bp, x_pre_flyby, x_post_flyby, ...
    state_pre_tcm1, epoch_tcm1, t_tcm2, epoch_flyby_ref, ...
    target_position, epoch_comet_flyby, c, dv_tcm);

% Convergence flags of the two ephemeris refinement stages.
% Used by main_organized to decide whether to plot the refined point.
out.success_pre  = info_pre.success;
out.success_post = info_post.success;

fprintf('\n--- Refinement convergence ---\n');
fprintf('  pre-flyby : %s (exitflag=%d, cv=%.2e)\n', ...
    choose_str(out.success_pre,  'OK', 'FAILED'), info_pre.exitflag,  info_pre.constrviolation);
fprintf('  post-flyby: %s (exitflag=%d, cv=%.2e)\n', ...
    choose_str(out.success_post, 'OK', 'FAILED'), info_post.exitflag, info_post.constrviolation);

fprintf('\nRefinement complete.\n');
end

%% Local helper
function s = choose_str(flag, s_true, s_false)
% choose_str  Return s_true if flag is true, otherwise s_false.
if flag, s = s_true; else, s = s_false; end
end







