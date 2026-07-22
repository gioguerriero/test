function [pareto_data, fig_pareto] = compute_pareto_front(global_results, c, S_halo, unstable_dir, eps_vel_ms)
% compute_pareto_front  Compute the true Pareto front (TOF vs total DV) over all
% solutions.
%
% For each solution it re-propagates the halo-to-moon and moon (back) arcs using
% the optimal parameters stored in s.x_opt, computes the effective DSM1 delta-v
% (velocity mismatch at the encounter), then builds the true Pareto front in the
% (total dv, total tof) plane using all solutions as individual points (no
% clustering by dv_max).
%
% Inputs:
%   global_results - cell array of solution structs from cr3bp_search
%   c              - constants struct (mu, Vstar, Tstar, ...)
%   S_halo         - halo orbit states array (for state_finder)
%   unstable_dir   - unstable directions for the initial perturbation
%   eps_vel_ms     - initial perturbation magnitude [m/s] (e.g. 15)
%
% Outputs:
%   pareto_data.sorted_solutions - all solutions sorted by TOF
%   pareto_data.dv_all           - effective total dv of each solution [km/s]
%   pareto_data.tof_all          - total tof of each solution [days]
%   pareto_data.dv_h2m_all       - h2m dv only (DSM1) [km/s]
%   pareto_data.dv_front         - DV points on the Pareto front [km/s]
%   pareto_data.tof_front        - TOF points on the Pareto front [days]
%   pareto_data.idx_front        - indices (in global_results) of the front points
%   fig_pareto                   - Pareto-front figure handle ([] if not created)

fig_pareto = [];   % figure handle (assigned later if the figure is created)

if isempty(global_results)
    warning('global_results is empty.');
    pareto_data.sorted_solutions = {};
    pareto_data.dv_all   = [];
    pareto_data.tof_all  = [];
    pareto_data.dv_front  = [];
    pareto_data.tof_front = [];
    pareto_data.idx_front = [];
    return;
end

n_sol = length(global_results);

%% Re-propagation: compute the effective dv_h2m for each solution
dv_injection_kms = eps_vel_ms * 1e-3;   % 15 m/s -> km/s
eps_vel_ad       = eps_vel_ms / (1e3 * c.Vstar);

dv_h2m_all = nan(1, n_sol);
dv_all     = nan(1, n_sol);
tof_all    = nan(1, n_sol);

opt = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);

fprintf('Re-propagating %d solutions to compute true dv_h2m...\n', n_sol);
for k = 1:n_sol
    s = global_results{k};

    % Optimal parameters
    moon_fpa          = s.x_opt(1);
    moon_out_of_plane = s.x_opt(2);
    halo_theta        = s.x_opt(3);
    tof_moon          = s.x_opt(4);
    tof_halo          = s.x_opt(5);

    moon_synodic   = s.match.moon_synodic;
    vinf_fm        = s.fmincon_vinf_results.vinf_fmincon;
    fpa_fm         = s.fmincon_vinf_results.fpa_fmincon;
    oop_fm         = s.fmincon_vinf_results.out_of_plane_fmincon;

    % v0_norm: reference v_inf magnitude (as in refinement_halo2moon)
    v0      = vinf_rotation(moon_synodic, vinf_fm, fpa_fm, oop_fm);
    v0_norm = norm(v0);

    try
        % --- Moon arc (back-propagation) ---
        v_inf = vinf_rotation(moon_synodic, v0_norm, moon_fpa, moon_out_of_plane);
        S0    = moon_synodic;
        S0(4:6) = S0(4:6) + v_inf ./ c.Vstar;
        [~, S_moon_mat] = ode45(@(t,S) CR3BP(t,S,c.mu), [0 -tof_moon], S0, opt);

        % --- Halo arc (forward-propagation) ---
        [initial_state, idx] = state_finder(halo_theta, S_halo);
        vu = unstable_dir(idx,:)';
        vu = vu / norm(vu);
        S0_pert = [initial_state(1:3); initial_state(4:6) + eps_vel_ad*vu];
        [~, S_halo_mat] = ode45(@(t,S) CR3BP(t,S,c.mu), [0 tof_halo], S0_pert, opt);

        % --- dv_h2m = velocity mismatch at rendezvous ---
        v_end_halo = S_halo_mat(end, 4:6);
        v_end_moon = S_moon_mat(end, 4:6);
        dv_h2m_kms = norm(v_end_halo - v_end_moon) * c.Vstar;

        dv_h2m_all(k) = dv_h2m_kms;
        dv_all(k)     = dv_injection_kms + dv_h2m_kms + s.dv2;
        tof_all(k)    = s.tof_total_days;
    catch ME
        warning('Solution %d: re-propagation failed (%s). Skipped.', k, ME.message);
    end
end

valid = ~isnan(dv_all);
n_valid = sum(valid);
fprintf('Re-propagation done: %d/%d solutions valid.\n', n_valid, n_sol);

%% Sort all solutions by total TOF (for the ranking table)
all_tof_sort = cellfun(@(s) s.tof_total_days, global_results);
[~, sort_idx] = sort(all_tof_sort, 'ascend');
pareto_data.sorted_solutions = global_results(sort_idx);

%% Print ranking table (now includes effective DV)
n_print = min(20, n_sol);
fprintf('\nAll solutions ranked by total TOF (%d shown of %d):\n', n_print, n_sol);
fprintf('  %-5s  %-12s  %-10s  %-10s  %-6s  %-10s  %-12s  %-12s\n', ...
    'Rank', 'TOF_tot[d]', 'TOF_h2m[d]', 'TOF_m2c[d]', 'q', 'vinf_ub', 'dv_h2m[km/s]', 'dv_tot[km/s]');
for ip = 1:n_print
    s = pareto_data.sorted_solutions{ip};
    k_orig = sort_idx(ip);
    fprintf('  %-5d  %-12.2f  %-10.2f  %-10.2f  %-6.2f  %-10.4f  %-12.4f  %-12.4f\n', ...
        ip, s.tof_total_days, s.tof_h2m_days, s.tof_m2c_days, s.q, s.vinf_ub, ...
        dv_h2m_all(k_orig), dv_all(k_orig));
end

%% Pareto front (non-dominance over all valid solutions)
idx_valid = find(valid);
dv_v   = dv_all(idx_valid);
tof_v  = tof_all(idx_valid);
n_v    = length(idx_valid);

is_nd = true(1, n_v);
for ip = 1:n_v
    for jp = 1:n_v
        if ip == jp, continue; end
        % jp dominates ip if it is <= on both objectives and < on at least one
        if dv_v(jp) <= dv_v(ip) && tof_v(jp) <= tof_v(ip) && ...
           (dv_v(jp) < dv_v(ip) || tof_v(jp) < tof_v(ip))
            is_nd(ip) = false;
            break;
        end
    end
end

idx_front_local = find(is_nd);
[dv_front, ord] = sort(dv_v(idx_front_local));
tof_tmp   = tof_v(idx_front_local);
tof_front = tof_tmp(ord);
idx_front = idx_valid(idx_front_local(ord));   % indices in the original global_results

pareto_data.dv_all      = dv_all;
pareto_data.tof_all     = tof_all;
pareto_data.dv_h2m_all  = dv_h2m_all;
pareto_data.dv_front    = dv_front;
pareto_data.tof_front   = tof_front;
pareto_data.idx_front   = idx_front;

fprintf('\nPareto front: %d non-dominated solutions.\n', length(dv_front));

%% Plot
if n_valid == 0, return; end

fig = figure('Name', 'Pareto Front — DV_total vs TOF', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.15 0.15 0.60 0.55]);
fig_pareto = fig;   % returned to the caller for saving as .fig
hold on; grid on; box on;

% All solutions (grey)
plot(tof_v, dv_v*1e3, 'o', ...
    'Color', [0.65 0.65 0.65], 'MarkerFaceColor', [0.65 0.65 0.65], ...
    'MarkerSize', 6, 'DisplayName', 'All solutions');

% Pareto front (blue), sorted by increasing TOF
[tof_front_plot, ord_plot] = sort(tof_front);
dv_front_plot = dv_front(ord_plot);
plot(tof_front_plot, dv_front_plot*1e3, 'b-o', 'LineWidth', 2, ...
    'MarkerFaceColor', 'b', 'MarkerSize', 9, 'DisplayName', 'Pareto front');

for ip = 1:length(dv_front_plot)
    text(tof_front_plot(ip) + 1.5, dv_front_plot(ip)*1e3, ...
        sprintf('%.1f m/s', dv_front_plot(ip)*1e3), 'FontSize', 9, 'Color', 'b');
end

xlabel('TOF_{total}  [days]',    'FontSize', 13);
ylabel('\DeltaV_{total}  [m/s]', 'FontSize', 13);
title('Pareto Front — \DeltaV_{total} vs TOF_{total}', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 11);
% NB: the figure is saved as .fig in the run folder by main_organized via the
%     'fig_pareto' handle. No file is written here.
end
