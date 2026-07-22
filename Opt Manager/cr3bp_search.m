function global_results = cr3bp_search(c, selected_comet, S_halo, unstable_dir, sp)
%CR3BP_SEARCH  CR3BP grid search + halo-to-moon refinement.
%
%  Inputs:
%    c              - constants struct
%    selected_comet - comet struct (must include comet_pos, epoch, theta,
%                     comet_pos_CR3BP, earth_pos)
%    S_halo         - halo orbit states matrix
%    unstable_dir   - unstable eigenvectors matrix
%    sp             - search_params struct (see main.m)
%
%  Output:
%    global_results - cell array of solution structs, one per saved solution

%% Unpack parameters
eps_vel_ms           = sp.eps_vel_ms;
h_min_earth          = sp.h_min_earth;
vinf_ub_vec          = sp.vinf_ub_vec;
q_vec                = sp.q_vec;
maximum_dv_vec       = sp.maximum_dv_vec;
N_save               = sp.N_save;
M                    = sp.M;
tmax_years_from_moon = sp.tmax_years_from_moon;
points_moon          = sp.points_moon;
limit_dist_km        = sp.limit_dist_km;
K_per_moon           = sp.K_per_moon;
div_tof_days         = sp.div_tof_days;
div_theta            = sp.div_theta;
k_dv_filter          = sp.k_dv_filter;
dv_floor_kms         = sp.dv_floor_kms;
max_matches          = sp.max_matches;
vinf.lower_bound     = sp.vinf_lower_bound;
min_days_between     = sp.min_days_between;   % [inj->DSM1, DSM1->flyby, flyby->DSM2, DSM2->comet] [days]

%% Build or load halo manifold tree
if sp.create_new_tree

    fprintf('Building halo manifold database (N=%d)...\n', sp.N);
    manifold_fn = build_halo_manifold_db(S_halo, unstable_dir, c, sp.N, ...
        eps_vel_ms, sp.tmax_years_from_halo, h_min_earth, sp.points_halo);

    fprintf('Building KD-tree...\n');
    Nh          = length(manifold_fn);
    all_halo_pos = cell2mat(arrayfun(@(s) s.state(:,1:3), manifold_fn(:), 'UniformOutput', false));
    all_halo_vel = cell2mat(arrayfun(@(s) s.state(:,4:6), manifold_fn(:), 'UniformOutput', false));
    n_pts_each   = arrayfun(@(s) size(s.state,1), manifold_fn(:));

    halo_tree_data.tree           = KDTreeSearcher(all_halo_pos);
    halo_tree_data.all_halo_vel   = all_halo_vel;
    halo_tree_data.halo_traj_id   = repelem((1:Nh)', n_pts_each);
    halo_tree_data.halo_time_id   = cell2mat(arrayfun(@(n) (1:n)', n_pts_each, 'UniformOutput', false));
    halo_tree_data.all_halo_theta = arrayfun(@(s) s.theta, manifold_fn(:));
    fprintf('KD-tree built: %d halo points.\n', size(all_halo_pos,1));

    save('halo_tree.mat',     'halo_tree_data');
    save('halo_manifold.mat', 'manifold_fn');

else
    load('halo_tree.mat');
    load('halo_manifold.mat');
end

% --- BEFORE the loops (right after loading/building the tree) ---
% S_halo_C       = parallel.pool.Constant(S_halo);
% unstable_dir_C = parallel.pool.Constant(unstable_dir);
% c_C            = parallel.pool.Constant(c);

% delete(gcp('nocreate'))

%% Search
n_vinf = length(vinf_ub_vec);
n_q    = length(q_vec);

fprintf('\n==========================================================\n');
fprintf('  GLOBAL SEARCH: %d vinf_ub x %d q = %d combinations\n', ...
    n_vinf, n_q, n_vinf * n_q);
fprintf('==========================================================\n\n');

global_results  = {};
global_counter  = 0;
min_tof_halo    = min_days_between(1) * 24*3600 / c.Tstar;   % min inj->DSM1 leg [adim]

for iv = 1:n_vinf
    for iq = 1:n_q

        vinf.upper_bound = vinf_ub_vec(iv);
        q                = q_vec(iq);

        fprintf('\n----------------------------------------------------------\n');
        fprintf('  [%d/%d] vinf_ub = %.2f km/s | q = %.2f\n', ...
            (iv-1)*n_q + iq, n_vinf*n_q, vinf.upper_bound, q);
        fprintf('----------------------------------------------------------\n');

        %% Moon-to-Comet optimization
        % Reuse an existing parallel pool instead of always calling parpool
        pool = gcp('nocreate');            % existing pool, or [] if none exists
        if isempty(pool)
            parpool('local', 8);           % create it only if none exists
        end

        try
            [result_ga_m2c, result_fmincon_struc_m2c, result_fmincon_modtof] = ...
                optimization_moon2comet(c, selected_comet.comet_pos, vinf, ...
                selected_comet.theta, selected_comet.comet_pos_CR3BP, ...
                selected_comet.epoch, q, min_days_between);
        catch ME
            fprintf('  optimization_moon2comet FAILED: %s\n', ME.message);
            continue;
        end
        delete(gcp('nocreate'))

        result_fmincon       = result_fmincon_modtof.bestfeasible.x;
        vinf_fmincon         = result_fmincon(1);
        theta_moon           = result_fmincon(2);
        fpa_fmincon          = result_fmincon(3);
        out_of_plane_fmincon = result_fmincon(4);
        dv2                  = result_fmincon_modtof.bestfeasible.fval;   % [km/s]
        tof_m2c              = result_fmincon(5);                          % [adim]

        %% Build Moon manifold
        manifold_back_fn = build_moon_manifold_db(c, theta_moon, vinf_fmincon, ...
            fpa_fmincon, out_of_plane_fmincon, M, tmax_years_from_moon, ...
            h_min_earth, points_moon);

        %% Matching
        matches_fn = find_manifold_matches(manifold_fn, manifold_back_fn, c, ...
            limit_dist_km, K_per_moon, div_tof_days, div_theta, ...
            halo_tree_data, min_tof_halo);

        if isempty(matches_fn)
            fprintf('  No matches found -> skip\n');
            continue;
        end
        fprintf('  Matches found: %d\n', length(matches_fn));

        matches_fn_orig = matches_fn;

        for i_dv = 1:length(maximum_dv_vec)

            matches_fn  = matches_fn_orig;
            maximum_dv  = maximum_dv_vec(i_dv);
            dv_remained = maximum_dv - eps_vel_ms*1e-3 - dv2;

            if dv_remained <= 0
                fprintf('  dv_remained = %.4f km/s -> skip (no budget)\n', dv_remained);
                continue;
            end
            fprintf('  dv2 = %.4f km/s | dv_remained = %.4f km/s\n', dv2, dv_remained);

            %% Filter matches by DV
            dv_norm_all = [matches_fn.dv_norm];
            dv_thresh   = max(k_dv_filter * dv_remained, dv_floor_kms);
            matches_fn  = matches_fn(dv_norm_all * c.Vstar <= dv_thresh);
            fprintf('  Matches after DV filter: %d\n', length(matches_fn));

            if isempty(matches_fn), continue; end

            %% Sort by total ToF ascending
            [~, ord] = sort([matches_fn.tof], 'ascend');
            matches_fn = matches_fn(ord);

            %% Greedy diversity filter
                
            relaxer = 1.2;

            tol_theta = relaxer * 15 * pi/180;               % [rad]
            tol_tof_h = relaxer * 20 * 24 * 3600 / c.Tstar; % [adim]
            tol_tof_m = relaxer * 20 * 24 * 3600 / c.Tstar; % [adim]
            tol_fpa   = relaxer * 20 * pi/180;               % [rad]
            tol_oop   = relaxer * 20 * pi/180;               % [rad]

            sel       = false(1, length(matches_fn));
            sel_theta = [];
            sel_tof_h = [];
            sel_tof_m = [];
            sel_fpa   = [];
            sel_oop   = [];

            for ic = 1:length(matches_fn)
                % if sum(sel) >= max_matches, break; end

                m        = matches_fn(ic);
                theta_ic = m.halo_theta;
                tof_h_ic = m.halo_time_of_encounter;
                tof_m_ic = abs(m.moon_time_of_encounter);
                fpa_ic   = m.moon_fpa;
                oop_ic   = m.moon_out_of_plane;

                is_diverse = true;
                for ks = 1:length(sel_theta)
                    if abs(theta_ic - sel_theta(ks)) < tol_theta && ...
                       abs(tof_h_ic - sel_tof_h(ks)) < tol_tof_h && ...
                       abs(tof_m_ic - sel_tof_m(ks)) < tol_tof_m && ...
                       abs(fpa_ic   - sel_fpa(ks))   < tol_fpa   && ...
                       abs(oop_ic   - sel_oop(ks))   < tol_oop
                        is_diverse = false;
                        break;
                    end
                end

                if is_diverse
                    sel(ic)   = true;
                    sel_theta = [sel_theta; theta_ic];
                    sel_tof_h = [sel_tof_h; tof_h_ic];
                    sel_tof_m = [sel_tof_m; tof_m_ic];
                    sel_fpa   = [sel_fpa;   fpa_ic];
                    sel_oop   = [sel_oop;   oop_ic];
                end
            end

            matches_fn = matches_fn(sel);
            fprintf('  Matches after diversity filter: %d\n', length(matches_fn));

            if isempty(matches_fn), continue; end

            %% Halo-to-Moon refinement
            fmincon_vinf_results.theta_moon           = theta_moon;
            fmincon_vinf_results.vinf_fmincon         = vinf_fmincon;
            fmincon_vinf_results.fpa_fmincon          = fpa_fmincon;
            fmincon_vinf_results.out_of_plane_fmincon = out_of_plane_fmincon;

            n_matches = length(matches_fn);
            matches_refined_tmp(n_matches) = struct( ...
                'x_opt', [], 'tof', [], 'exitflag', [], 'guess_index', []);
            tof_refined_vec = nan(1, n_matches);

            for i = 1:n_matches % (parfor could be used here)
                try
                    [result_i, tof_i, exitflag_i] = refinement_halo2moon(c, matches_fn(i), ...
                        dv_remained, vinf_fmincon, fpa_fmincon, out_of_plane_fmincon, ...
                        h_min_earth, eps_vel_ms, S_halo, unstable_dir, fmincon_vinf_results, min_days_between);

                    % the variables needed for the parpool are set up above
                    % [result_i, tof_i, exitflag_i] = refinement_halo2moon(c_C.Value, matches_fn(i), ...
                    %     dv_remained, vinf_fmincon, fpa_fmincon, out_of_plane_fmincon, ...
                    %     h_min_earth, eps_vel_ms, S_halo_C.Value, unstable_dir_C.Value, ...
                    %     fmincon_vinf_results, min_days_between);

                    if exitflag_i > 0
                        matches_refined_tmp(i).x_opt      = result_i;
                        matches_refined_tmp(i).tof        = tof_i;
                        matches_refined_tmp(i).exitflag   = exitflag_i;
                        matches_refined_tmp(i).guess_index= i;
                        tof_refined_vec(i)                = tof_i;
                    end
                catch
                end
            end

            valid_idx       = ~isnan(tof_refined_vec);
            matches_refined = matches_refined_tmp(valid_idx);
            tof_refined_vec = tof_refined_vec(valid_idx);
            clear matches_refined_tmp;

            fprintf('  Converged: %d / %d\n', length(matches_refined), n_matches);
            if isempty(matches_refined), continue; end

            %% Save top N_save solutions by total TOF
            tof_total_vec = tof_refined_vec + tof_m2c;
            [~, sort_ord] = sort(tof_total_vec, 'ascend');
            mr_sorted = matches_refined(sort_ord);
            n_to_save = min(N_save, length(mr_sorted));


            for is = 1:n_to_save
                global_counter = global_counter + 1;

                sol.x_opt                 = mr_sorted(is).x_opt;
                sol.tof_h2m               = mr_sorted(is).tof;
                sol.tof_m2c               = tof_m2c;
                sol.tof_total             = mr_sorted(is).tof + tof_m2c;
                sol.tof_total_days        = (mr_sorted(is).tof + tof_m2c) * c.Tstar / 86400;
                sol.tof_h2m_days          = mr_sorted(is).tof * c.Tstar / 86400;
                sol.tof_m2c_days          = tof_m2c * c.Tstar / 86400;
                sol.guess_index           = mr_sorted(is).guess_index;
                sol.vinf_ub               = vinf.upper_bound;
                sol.q                     = q;
                sol.dv2                   = dv2;
                sol.dv_remained           = dv_remained;
                sol.maximum_dv            = maximum_dv;
                sol.match                 = matches_fn(mr_sorted(is).guess_index);
                % sol.result_ga_m2c         = result_ga_m2c;
                sol.result_fmincon        = result_fmincon;
                % sol.result_fmincon_modtof = result_fmincon_modtof;
                sol.fmincon_vinf_results  = fmincon_vinf_results;

                global_results{global_counter} = sol;
            end

            fprintf('  Top %d saved (best TOF = %.2f days)\n', n_to_save, ...
                tof_total_vec(sort_ord(1)) * c.Tstar / 86400);
        end
    end
end

fprintf('\n  Total solutions stored: %d\n', global_counter);
end
