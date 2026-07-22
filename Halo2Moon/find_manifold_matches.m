function matches = find_manifold_matches(manifold, manifold_back, c, limit_dist_km, K_per_moon, div_tof_days, div_theta, halo_tree_data, min_tof_halo)
% find_manifold_matches  Match Halo and Moon manifold points using a KD-tree,
% rangesearch, and a greedy diversity filter.
%
% For each Moon trajectory:
%   1. rangesearch: find all Halo points within limit_dist
%   2. Sort by increasing delta-V
%   3. Greedy selection: keep up to K_per_moon mutually diverse matches
%      (diverse = |dtof_moon| > div_tof OR |dtof_halo| > div_tof OR |dtheta| > div_theta)
%
% Inputs:
%   manifold       - Halo manifold struct array
%   manifold_back  - Moon (back-propagated) manifold struct array
%   c              - constants struct
%   limit_dist_km  - distance threshold [km]
%   K_per_moon     - max matches per Moon trajectory (default 5)
%   div_tof_days   - min tof difference for diversity [days] (default 20)
%   div_theta      - min theta difference for diversity [rad] (default 0.3)
%   halo_tree_data - (optional) precomputed Halo KD-tree struct; built internally
%                    from manifold if omitted. Fields: tree, all_halo_vel,
%                    halo_traj_id, halo_time_id, all_halo_theta
%   min_tof_halo   - minimum Halo time of flight to accept a match [non-dim]
%
% Outputs:
%   matches        - struct array of matches, ranked by increasing delta-V

% ---------------------- DEFAULT PARAMETERS ---------------------- %
if nargin < 5 || isempty(limit_dist_km),  limit_dist_km = 30000;  end
if nargin < 6 || isempty(K_per_moon),     K_per_moon = 5;         end
if nargin < 7 || isempty(div_tof_days),   div_tof_days = 20;      end
if nargin < 8 || isempty(div_theta),      div_theta = 0.3;        end

limit_dist = limit_dist_km / c.Lstar;
div_tof    = div_tof_days * 24 * 3600 / c.Tstar;  % non-dimensional

Nh = length(manifold);
Nm = length(manifold_back);

% ---------------------- GLOBAL KD-TREE (all Halo points) --------- %
if nargin >= 8 && ~isempty(halo_tree_data)
    % Use the precomputed tree passed in from outside
    tree           = halo_tree_data.tree;
    all_halo_vel   = halo_tree_data.all_halo_vel;
    halo_traj_id   = halo_tree_data.halo_traj_id;
    halo_time_id   = halo_tree_data.halo_time_id;
    all_halo_theta = halo_tree_data.all_halo_theta;
    fprintf('Using pre-built KD-tree.\n');
else
    % Build the tree internally
    all_halo_pos = cell2mat(arrayfun(@(s) s.state(:,1:3), manifold(:), 'UniformOutput', false));
    all_halo_vel = cell2mat(arrayfun(@(s) s.state(:,4:6), manifold(:), 'UniformOutput', false));
    n_pts_each   = arrayfun(@(s) size(s.state,1), manifold(:));
    halo_traj_id = repelem((1:Nh)', n_pts_each);
    halo_time_id = cell2mat(arrayfun(@(n) (1:n)', n_pts_each, 'UniformOutput', false));
    all_halo_theta = arrayfun(@(s) s.theta, manifold(:));

    fprintf('Building global KD-tree (%d Halo points)...\n', size(all_halo_pos,1));
    tree = KDTreeSearcher(all_halo_pos);
end

% ---------------------- PRE-EXTRACT MOON DATA -------------------- %
moon_pos = cell(Nm,1);
moon_vel = cell(Nm,1);
for j = 1:Nm
    moon_pos{j} = manifold_back(j).state(:,1:3);
    moon_vel{j} = manifold_back(j).state(:,4:6);
end

% ---------------------- MAIN MATCHING LOOP ----------------------- %
results_cell = {};
counter = 0;

for j = 1:Nm
    if mod(j, 50) == 0 || j == Nm
        fprintf('Moon traj %d / %d\n', j, Nm);
    end

    pos_m = moon_pos{j};
    vel_m = moon_vel{j};
    t_moon_j = manifold_back(j).t;

    % Find all Halo points within limit_dist for each Moon point
    [idx_all, dist_all] = rangesearch(tree, pos_m, limit_dist);

    % ---- Vectorized collection of all candidates ---- %
    n_neighbors = cellfun(@numel, idx_all);   % number of matches per Moon point
    total_cand  = sum(n_neighbors);

    if total_cand == 0, continue; end

    % rangesearch returns row vectors per cell -> horzcat + transpose to get a column
    cand_halo_gid = [idx_all{:}]';                                  % [total_cand x 1]
    cand_dist     = [dist_all{:}]';                                 % [total_cand x 1]
    cand_moon_idx = repelem((1:size(pos_m,1))', n_neighbors);       % [total_cand x 1]

    dv_vecs = all_halo_vel(cand_halo_gid, :) - vel_m(cand_moon_idx, :);
    cand_dv = vecnorm(dv_vecs, 2, 2);                              % [total_cand x 1]

    % ---- Sort by increasing delta-V ---- %
    [~, sort_order] = sort(cand_dv, 'ascend');
    cand_moon_idx = cand_moon_idx(sort_order);
    cand_halo_gid = cand_halo_gid(sort_order);
    cand_dv       = cand_dv(sort_order);
    cand_dist     = cand_dist(sort_order);

    % ---- Greedy selection with a diversity criterion ---- %
    % Keep the first (min dv), then for each next candidate require it to be
    % different from ALL already-kept ones:
    %   |dtof_moon| > div_tof  OR  |dtof_halo| > div_tof  OR  |dtheta| > div_theta

    sel_idx     = [];   % indices into the sorted candidates
    sel_theta   = [];
    sel_tof_h   = [];
    sel_tof_m   = [];

    for ic = 1:length(cand_dv)
        if length(sel_idx) >= K_per_moon, break; end

        gid   = cand_halo_gid(ic);
        h_idx = halo_traj_id(gid);
        h_tid = halo_time_id(gid);
        m_idx = cand_moon_idx(ic);

        theta_ic    = all_halo_theta(h_idx);
        tof_halo_ic = manifold(h_idx).t(h_tid);
        if tof_halo_ic < min_tof_halo, continue; end  % skip points too close to the start
        tof_moon_ic = t_moon_j(m_idx);

        % Check diversity against all already-selected matches
        is_diverse = true;
        for is = 1:length(sel_idx)
            d_theta = abs(theta_ic    - sel_theta(is));
            d_tof_h = abs(tof_halo_ic - sel_tof_h(is));
            d_tof_m = abs(tof_moon_ic - sel_tof_m(is));

            % If NONE of the three exceeds its threshold -> too similar
            if d_theta < div_theta && d_tof_h < div_tof && d_tof_m < div_tof
                is_diverse = false;
                break;
            end
        end

        if is_diverse
            sel_idx   = [sel_idx;   ic];
            sel_theta = [sel_theta; theta_ic];
            sel_tof_h = [sel_tof_h; tof_halo_ic];
            sel_tof_m = [sel_tof_m; tof_moon_ic];
        end
    end

    % ---- Save the selected matches ---- %
    for is = 1:length(sel_idx)
        ic    = sel_idx(is);
        gid   = cand_halo_gid(ic);
        h_idx = halo_traj_id(gid);
        h_tid = halo_time_id(gid);
        m_idx = cand_moon_idx(ic);

        counter = counter + 1;

        t_moon_full = manifold_back(j).t;
        t_halo_full = manifold(h_idx).t;

        result = struct();
        result.moon_index             = j;
        result.halo_index             = h_idx;
        result.moon_fpa               = manifold_back(j).fpa;
        result.moon_out_of_plane      = manifold_back(j).out_of_plane;
        result.halo_theta             = manifold(h_idx).theta;
        result.moon_time_of_encounter = t_moon_full(m_idx);
        result.halo_time_of_encounter = t_halo_full(h_tid);
        result.dv_norm                = cand_dv(ic);
        result.tof                    = abs(t_moon_full(m_idx)) + t_halo_full(h_tid);
        result.min_dist               = cand_dist(ic);
        result.halo_state_point       = manifold(h_idx).state(h_tid,:);
        result.moon_state_point       = manifold_back(j).state(m_idx,:);
        result.halo_traj              = manifold(h_idx);
        result.moon_traj              = manifold_back(j);
        result.moon_synodic           = manifold_back(j).moon_synodic;

        results_cell{counter} = result;
    end
end

% ---------------------- FINAL OUTPUT ----------------------------- %
if counter == 0
    disp("No match found.");
    matches = struct([]);
else
    matches = [results_cell{:}];

    % Rank by delta-V ascending
    dv_list = [matches.dv_norm]';
    [~, order] = sort(dv_list, 'ascend');
    matches = matches(order);

    fprintf('Matching completed. %d matches saved (K=%d per Moon traj).\n', counter, K_per_moon);
end

end
