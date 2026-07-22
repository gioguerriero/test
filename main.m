%% ================================================================
%  MAIN — Comet Mission Trajectory Optimization
%
%  Pipeline:
%    1. CR3BP global search  ->  global_results
%    2. Pareto front analysis
%    3. Solution selection (Pareto front)
%    4. Ephemeris refinement  ->  out
%    5. Final trajectory plot
% ================================================================
clear; clc; close all;

%% Paths
project_root = fileparts(mfilename('fullpath'));  % folder containing this main script
cd(project_root);

addpath('Auxiliar/');
addpath('Moon2Comet/');
addpath('Halo2Moon/');
addpath('CR3BP 2 Ephemeris Sun-Earth/');
addpath('Full_ephemeris_conversion/');
addpath('Full_ephemeris_conversion/Refinement post flyby/');
addpath('Full_ephemeris_conversion/Refinement pre flyby/');
addpath('Full_ephemeris_conversion/Refinement global post flyby/');
addpath('Full_ephemeris_conversion/tcm2tcm/');
addpath('Opt Manager/');
addpath('Results/');

startup
cspice_furnsh({'kernels/sckernel.tm'});

%% Physical constants
c.G      = 6.67430e-20;           % [km^3/kg/s^2]
c.Lstar  = 1.496e+8;              % [km]
c.Tstar  = 5021870.4424055;       % [s]
c.Vstar  = 29.806;                % [km/s]
c.mSun   = getAstroConstants('Sun',   'mass');
c.mEarth = getAstroConstants('Earth', 'mass');
c.mMoon  = 7.34767309e22;
c.rSun   = 695700;                % [km]
c.rEarth = 6378;                  % [km]
c.rMoon  = 1737;                  % [km]
c.rMoon_ad = 384400 / c.Lstar;
c.mu     = c.mEarth / (c.mSun + c.mEarth);
c.muMoon = 4902.800066;           % [km^3/s^2]

%% ================================================================
%  USER-CONFIGURABLE PARAMETERS
% ================================================================

% ---- CR3BP sweep ---------------------------------------------
vinf_ub_vec    = [0.7];   % vinf upper bound(s) [km/s]
q_vec          = [0.442857142857143];    % q parameter for optimization_moon2comet
maximum_dv_vec = [0.85];    % total DV budget(s) [km/s]                               
N_save         = 20;        % top solutions saved per (vinf_ub, q, dv_max)

% vinf_ub_vec    = linspace(0.7, 0.85, 5);   % vinf upper bound(s) [km/s]
% q_vec          = linspace(0.1, 0.45, 10);    % q parameter for optimization_moon2comet
% maximum_dv_vec = linspace(0.1, 0.85, 6);    % total DV budget(s) [km/s]

% ---- Moon manifold -------------------------------------------
M                    = 2000;
tmax_years_from_moon = 1;
points_moon          = 5000;

% ---- Matching ------------------------------------------------
limit_dist_km = 20000;
K_per_moon    = 10;
div_tof_days  = 30;
div_theta     = deg2rad(30);

% ---- Match filtering -----------------------------------------
k_dv_filter  = 2;     % keep matches with dv < max(k * dv_remained, dv_floor)
dv_floor_kms = 0.4;   % [km/s] minimum DV threshold
max_matches  = 500;   % max matches passed to refinement

% ---- Ephemeris refinement ------------------------------------
max_dv            = 850;    % total DV budget [m/s]
hp                = 750;    % min lunar flyby altitude [km]
max_dv_inj        = 150;    % max injection DV [m/s]
multiple_shooting = 1;      % 1 = multiple shooting, 0 = single shooting
use_fast_refinement = 0;    % 0 = slow, 1 = fast

% Multiple shooting: one node every N days, for each of the 4 phases
%   [Halo->DSM1, DSM1->Flyby, Flyby->DSM2, DSM2->Comet]
days_per_node = [5 5 20 60];  % [days/node] per phase (to be tuned!)


% ---- Maneuver spacing (minimum days between events) ----------
% Minimum spacing [days] between consecutive events. Applied to BOTH the
% CR3BP global search AND the ephemeris refinement (optimization_ephe).
%   order: [inj->DSM1, DSM1->flyby, flyby->DSM2, DSM2->comet]
min_days_between = [3.5 3.5 3.5 3.5];   % [days]

%% ================================================================
%  COMET DEFINITIONS
%  Target catalogue. Each comet holds its heliocentric ECLIPJ2000 encounter
%  position [km] and epoch [ET s]. One entry is selected below and enriched
%  with the Earth state to build the CR3BP encounter geometry.
% ================================================================
C2023X1.comet_pos = [-1.1791e+8, 1.0825e+8, -2.1045e+5];
C2023X1.epoch     = 753263052;
C2023X1.name      = 'C-2023 X1'; 

C2013US10.comet_pos = [-125534119.984888, -13360764.5867413, 0.00683142548405158];
C2013US10.epoch     = 5807.96341905708 * 86400;
C2013US10.name      = 'C-2013 US10';

C2001Q4.comet_pos = [-123707083.407748, -72853679.6310874, 0.00397940697075683];
C2001Q4.epoch     = 1595.65383997182 * 86400;
C2001Q4.name      = 'C-2001 Q4';

C2008A1.comet_pos = [22000353.1090486, -159925043.720566, 0.000403648958354097];
C2008A1.epoch     = 3202.88049840330 * 86400;
C2008A1.name      = 'C-2008 A1';

C2013K1.comet_pos = [83618036.6060367, -115939662.210283, -9.94606334643322e-05];
C2013K1.epoch     = 4903.44378260833 * 86400;
C2013K1.name      = 'C-2013 K1';

C2015P3.comet_pos = [-3605089.65918768, -127156821.467427, 0.00241260529037390];
C2015P3.epoch     = 5709.36441765484 * 86400;
C2015P3.name      = 'C-2015 P3';

C2016VZ18.comet_pos = [-23038534.2141767, 149609028.629560, 0.00320507358901523];
C2016VZ18.epoch     = 6250.56247394526 * 86400;
C2016VZ18.name      = 'C-2016 VZ18'; 

C2024Y1.comet_pos = [60687786.9561721	-149785951.733452	-0.00240193589752380];
C2024Y1.epoch     = 9058.62066774852 * 86400;
C2024Y1.name      = 'C-2024 Y1'; 

C2019U6.comet_pos = [-83524908.6312209	-120256148.191715	-0.000103627814496576];
C2019U6.epoch     = 7494.29178284879 * 86400;
C2019U6.name      = 'C-2019 U6';

Synthetic1.comet_pos = [-61872301.7008814	-139144956.367171	-9.92433567673097e-07];
Synthetic1.epoch = 1.060917897278888e+09;
Synthetic1.name = 'Synthetic1';

Synthetic2.comet_pos = [-2.347967943883641e+07,1.490052698995286e+08,7.488288527911041e-07];
Synthetic2.epoch = 1.013411179613203e+09;
Synthetic2.name = 'Synthetic2';

Synthetic3.comet_pos = [-9.800464865785849e+07,1.535075228429915e+08,-5.774306732388979e-08];
Synthetic3.epoch = 1.045040734509147e+09; 
Synthetic3.name = 'Synthetic3';

% 3I/ATLAS 
C2025N1.comet_pos = [-2.0851e+8 -9.2347e+05 1.094e+07];
C2025N1.earth_pos = [1.0505e+08 1.0463e+08 0];
C2025N1.epoch = 9441.5 * 86400; 
C2025N1.name = "3IAtlas";

% Unusual comet intercepted out of the ecliptic plane
C2025R2.comet_pos = [147155133.376193 36074696.9451916 7308657.60371978];
C2025R2.epoch = 9425.7 * 86400;
C2025R2.name = 'C-2025 R2';

% ---- Select comet and enrich struct --------------------------
selected_comet = Synthetic1;   % comet processed by the pipeline

% Earth state at the encounter epoch defines the synodic frame orientation
selected_comet.earth_state = cspice_spkezr('EARTH', selected_comet.epoch, 'ECLIPJ2000', 'NONE', 'SUN');
selected_comet.earth_pos   = selected_comet.earth_state(1:3)';
selected_comet.theta       = atan2(selected_comet.earth_pos(2), selected_comet.earth_pos(1));

% Rotate the comet position into the CR3BP frame and normalise by the
% instantaneous Sun-Earth distance
Rz = [ cos(selected_comet.theta)  sin(selected_comet.theta)  0;
      -sin(selected_comet.theta)  cos(selected_comet.theta)  0;
       0                           0                           1];
comet_CR3BP                    = Rz * selected_comet.comet_pos';
comet_enc_Lstar                = norm(selected_comet.earth_pos(1:2));
selected_comet.comet_pos_CR3BP = comet_CR3BP ./ comet_enc_Lstar;

safe_name = strrep(strrep(selected_comet.name, '-', ''), ' ', '_');

%% ================================================================
%  RESULTS FOLDER — a new sub-folder is created for every run
%  Results/<CometName>_runX   (X = first free integer)
% ================================================================
results_root = fullfile(project_root, 'Results');
if ~exist(results_root, 'dir'); mkdir(results_root); end

run_prefix    = sprintf('%s_run', safe_name);
existing_runs = dir(fullfile(results_root, [run_prefix '*']));
used_nums     = [];
for ir = 1:numel(existing_runs)
    if existing_runs(ir).isdir
        n = sscanf(existing_runs(ir).name(numel(run_prefix)+1:end), '%d');
        if ~isempty(n); used_nums(end+1) = n(1); end %#ok<AGROW>
    end
end
if isempty(used_nums), run_num = 1; else, run_num = max(used_nums) + 1; end
run_dir = fullfile(results_root, sprintf('%s%d', run_prefix, run_num));
mkdir(run_dir);
fprintf('Results for this run will be saved to: %s\n', run_dir);

%% Load Halo data (used by both search and refinement)
load('S_halo.mat');
load('unstable_dir.mat');

%% Plot the comet positions in the synodic plane
comets = {C2023X1, C2013K1, C2019U6, C2024Y1, Synthetic1, Synthetic2, Synthetic3};  % comets to plot in the non-dimensional synodic frame

plot_comets_synodic(comets, c.mu, 'PlotManifolds', true, ...
    'S_halo', S_halo, 'UnstableDir', unstable_dir, 'C', c);

%% ================================================================
%  CR3BP GLOBAL SEARCH
%  Sweep the (vinf_ub, q, dv_max) grid to find low-cost Halo->Moon->Comet
%  transfers in the CR3BP. All search settings are gathered into
%  search_params and passed to cr3bp_search, which returns the candidate
%  solutions in global_results.
% ================================================================
search_params.vinf_ub_vec          = vinf_ub_vec;
search_params.q_vec                = q_vec;
search_params.maximum_dv_vec       = maximum_dv_vec;
search_params.N_save               = N_save;
search_params.vinf_lower_bound     = 0.2;       % [km/s]
search_params.eps_vel_ms           = 15;        % injection DV [m/s]
search_params.h_min_earth          = 10000;     % [km]
search_params.create_new_tree      = 0;         % 1 to rebuild halo tree
search_params.N                    = 3500;      % halo manifold points
search_params.tmax_years_from_halo = 1;
search_params.points_halo          = 5000;
search_params.M                    = M;
search_params.tmax_years_from_moon = tmax_years_from_moon;
search_params.points_moon          = points_moon;
search_params.limit_dist_km        = limit_dist_km;
search_params.K_per_moon           = K_per_moon;
search_params.div_tof_days         = div_tof_days;
search_params.div_theta            = div_theta;
search_params.k_dv_filter          = k_dv_filter;
search_params.dv_floor_kms         = dv_floor_kms;
search_params.max_matches          = max_matches;
search_params.min_days_between     = min_days_between;

tic;
global_results = cr3bp_search(c, selected_comet, S_halo, unstable_dir, search_params);
fprintf('CR3BP search completed in %.1f min\n', toc/60);

if isempty(global_results)
    error('No solutions found in CR3BP search. Adjust sweep parameters.');
end

% ---- Save this run's results: search_params + global_results ----
save(fullfile(run_dir, 'search_params.mat'),  'search_params');
save(fullfile(run_dir, 'global_results.mat'), 'global_results');

%% ================================================================
%  PARETO FRONT
%  Build the TOF vs total-DV trade-off across all solutions. Inspect the
%  printed table to choose which front solution(s) to refine in the loop
%  below. The result is also repackaged (sorted by TOF) into the
%  pareto_front struct used below.
% ================================================================
[pareto_data, fig_pareto] = compute_pareto_front(global_results, c, S_halo, unstable_dir, search_params.eps_vel_ms);

% ---- Save pareto_data + the Pareto-front figure as .fig ----
save(fullfile(run_dir, 'pareto_data.mat'), 'pareto_data');
if ~isempty(fig_pareto) && isgraphics(fig_pareto)
    savefig(fig_pareto, fullfile(run_dir, 'pareto_front.fig'));
end

% ---- Repackage only the Pareto-front solutions, sorted by TOF ----
% pareto_data.idx_front / tof_front / dv_front are ordered by increasing dV;
% here everything is re-sorted by increasing TOF. pareto_front collects:
%   .sols     - cell array with the full data of each front solution
%   .idx      - corresponding indices in global_results
%   .tof_days - vector of all TOFs [days]
%   .dv_kms   - vector of all total dVs [km/s]
%   .dv_ms    - vector of all total dVs [m/s]
[tof_sorted, ord] = sort(pareto_data.tof_front(:), 'ascend');

idx_sorted    = pareto_data.idx_front(ord);
dv_sorted_kms = pareto_data.dv_front(ord);

pareto_front.idx      = idx_sorted(:);          % column
pareto_front.tof_days = tof_sorted(:);          % column [days]
pareto_front.dv_kms   = dv_sorted_kms(:);       % column [km/s]
pareto_front.dv_ms    = dv_sorted_kms(:) * 1e3; % column [m/s]
pareto_front.sols     = global_results(pareto_front.idx);
pareto_front.sols     = pareto_front.sols(:);   % column of cells

% ---- Save the pareto_front struct ----
save(fullfile(run_dir, 'pareto_front.mat'), 'pareto_front');

%% Plot a chosen Pareto-front solution (CR3BP)
% Visualise the CR3BP trajectory of one Pareto-front solution, selected by its
% rank along the front. The front is sorted by increasing TOF, so
% pareto_rank = 1 is the fastest (lowest-TOF) solution, 2 the next, and so on.
pareto_rank = 1;   % <-- choose which Pareto-front solution to plot (1 = lowest TOF)

if pareto_rank < 1 || pareto_rank > numel(pareto_front.idx)
    warning(['pareto_rank = %d is out of range: the Pareto front has %d ' ...
             'solutions. Skipping the plot.'], pareto_rank, numel(pareto_front.idx));
else
    k_global = pareto_front.idx(pareto_rank);   % corresponding index in global_results
    plot_paper_trajectory(global_results, k_global, S_halo, c, selected_comet, ...
        unstable_dir, search_params.eps_vel_ms);
    fprintf('Plotted Pareto-front solution rank %d/%d (global idx %d): TOF=%.2f d | DV=%.1f m/s\n', ...
        pareto_rank, numel(pareto_front.idx), k_global, ...
        pareto_front.tof_days(pareto_rank), pareto_front.dv_ms(pareto_rank));
end



%% ================================================================
%  EPHEMERIS REFINEMENT — loop over all Pareto-front solutions
%  Each CR3BP solution is refined in the full ephemeris model. Node counts
%  (k_vec) are set per phase from the segment durations, then run_refinement
%  produces the converged trajectory (out) for every front solution.
% ================================================================
refine_params.max_dv            = max_dv;
refine_params.hp                = hp;
refine_params.max_dv_inj        = max_dv_inj;
refine_params.multiple_shooting = multiple_shooting;
refine_params.eps_vel_ms        = search_params.eps_vel_ms;
refine_params.t_tcm1            = 5;    % [days] before flyby
refine_params.use_fast          = use_fast_refinement;
refine_params.min_days_between  = min_days_between;

% ---- Save the refinement parameters ----
save(fullfile(run_dir, 'refine_params.mat'), 'refine_params');

N_front = numel(pareto_front.idx);

refined_front(N_front) = struct( ...
    'rank',           [], ...
    'idx',            [], ...
    'out',            [], ...
    'tof_total_days', [], ...
    'dv_total_ms',    [], ...
    'success',        [], ...
    'err_msg',        '');

for ii = 1:N_front
%for ii = [1] % if specific solution refinement is preferred

    k_orig = pareto_front.idx(ii);           % index in the original global_results
    best_i = pareto_front.sols{ii};          % solution (sorted by increasing TOF)

    % ---- Dynamic k_vec: one node every days_per_node days, per phase ----
    % Duration of the 4 segments [days] derived from the solution parameters.
    Tstar_d  = c.Tstar / 86400;              % days per non-dimensional time unit
    seg_days = [ best_i.x_opt(5); ...                                        % Halo  -> DSM1
                 best_i.x_opt(4); ...                                        % DSM1  -> Flyby
                 best_i.result_fmincon(5) *  best_i.result_fmincon(6); ...   % Flyby -> DSM2
                 best_i.result_fmincon(5) * (1 - best_i.result_fmincon(6)) ] * Tstar_d;
    k_vec_i  = max(1, round(seg_days(:).' ./ days_per_node));  % integers >= 1
    refine_params.k_vec = k_vec_i;

    fprintf('\n################################################################\n');
    fprintf('  Pareto solution %d/%d (global idx %d): TOF=%.2f d | vinf_ub=%.2f | q=%.2f | dv_max=%.0f m/s\n', ...
        ii, N_front, k_orig, best_i.tof_total_days, best_i.vinf_ub, best_i.q, best_i.maximum_dv*1e3);
    fprintf('################################################################\n');
    fprintf('  segments [d]: [%.1f %.1f %.1f %.1f] -> k_vec = [%d %d %d %d]\n', ...
        seg_days(1), seg_days(2), seg_days(3), seg_days(4), ...
        k_vec_i(1), k_vec_i(2), k_vec_i(3), k_vec_i(4));

    refined_front(ii).rank = k_orig;
    refined_front(ii).idx  = ii;

    tic;
    try
        out_i = run_refinement(c, best_i, selected_comet, S_halo, unstable_dir, refine_params);

        tof_total_days_i = (out_i.comet.epoch - out_i.departure.epoch) / 86400;
        dv_total_ms_i    = out_i.injection.norm_ms + out_i.dsm1.dv.norm_ms ...
                         + out_i.tcm1.dv.norm_ms   + out_i.tcm2.dv.norm_ms ...
                         + out_i.dsm2.dv.norm_ms;

        refined_front(ii).out            = out_i;
        refined_front(ii).tof_total_days = tof_total_days_i;
        refined_front(ii).dv_total_ms    = dv_total_ms_i;
        refined_front(ii).success        = true;
        refined_front(ii).err_msg        = '';

        fprintf('  -> refined: TOF=%.2f d | DV_tot=%.2f m/s  (%.1f min)\n', ...
            tof_total_days_i, dv_total_ms_i, toc/60);
    catch ME
        refined_front(ii).success = false;
        refined_front(ii).err_msg = ME.message;
        fprintf('  !! refinement failed: %s  (%.1f min)\n', ME.message, toc/60);
    end
end

% ---- Save the refined Pareto front ----
save(fullfile(run_dir, 'refined_front.mat'), 'refined_front');
fprintf('\nRun results saved to: %s\n', run_dir);

%% Close parallel pool
delete(gcp('nocreate'));

%% ================================================================
%  PLOT — CR3BP vs Refined Pareto front (TOF-DV plane, m/s)
%  Overlay the CR3BP front and the refined (ephemeris) front, linking each
%  CR3BP solution to its refined counterpart to show the cost/time shift.
% ================================================================
tof_cr3bp = pareto_front.tof_days(:);            % [days]  (sorted by TOF)
dv_cr3bp  = pareto_front.dv_ms(:);               % [m/s]   (sorted by TOF)

% Keep only solutions where refinement did not error AND both the pre-flyby
% and post-flyby optimizations converged.
ok = arrayfun(@(s) ~isempty(s.success) && s.success && ...
    ~isempty(s.out) && s.out.success_pre && s.out.success_post, refined_front);
tof_ref = arrayfun(@(s) s.tof_total_days, refined_front(ok));
dv_ref  = arrayfun(@(s) s.dv_total_ms,    refined_front(ok));

figure('Name', 'Pareto front: CR3BP vs Refined', 'Color', 'w');
hold on; grid on; box on;

% Segments linking each pair (CR3BP -> refined)
idx_ok = find(ok);
h_link = [];
for jj = 1:numel(idx_ok)
    ii = idx_ok(jj);
    h_link = plot([tof_cr3bp(ii), refined_front(ii).tof_total_days], ...
                  [dv_cr3bp(ii),  refined_front(ii).dv_total_ms], ...
                  '-', 'Color', [0.5 0.5 0.5 0.6], 'LineWidth', 0.8);
end

% Pareto front CR3BP
h_cr3bp = plot(tof_cr3bp, dv_cr3bp, 'o-', ...
    'Color', [0 0.4470 0.7410], 'MarkerFaceColor', [0 0.4470 0.7410], ...
    'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'CR3BP');

% Pareto front Refined
h_ref = plot(tof_ref, dv_ref, 's-', ...
    'Color', [0.8500 0.3250 0.0980], 'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
    'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Refined (ephemeris)');

xlabel('Total TOF [days]');
ylabel('Total \DeltaV [m/s]');
title(sprintf('Pareto front: CR3BP vs Refined — %s', selected_comet.name));
if ~isempty(h_link)
    legend([h_cr3bp, h_ref, h_link], {'CR3BP', 'Refined (ephemeris)', 'CR3BP \rightarrow Refined'}, ...
        'Location', 'best');
else
    legend([h_cr3bp, h_ref], 'Location', 'best');
end
hold off;

%% Plot the user-selected solution
% Choose which refined (full-ephemeris) solution to plot by changing k.
% k is the index in refined_front (same 'ii' used in the refinement loop,
% i.e. the rank along the Pareto front sorted by TOF).
k = 1;   % <-- change the index of the solution to visualise here

if k < 1 || k > numel(refined_front)
    error('k = %d out of range: refined_front has %d elements.', k, numel(refined_front));
elseif isempty(refined_front(k).out)
    error(['Solution k = %d was not refined ' ...
           '(refined_front(%d).out is empty). Choose a k included in the loop.'], k, k);
else
    out_k = refined_front(k).out;

    if isempty(refined_front(k).success) || ~refined_front(k).success
        warning('Solution k = %d did not converge: plotting it anyway.', k);
    end

    plot_full_trajectory(out_k, S_halo, c, selected_comet.name);

    fprintf('Plotted solution k=%d (global idx %d): TOF=%.2f d | DV_tot=%.2f m/s\n', ...
        k, refined_front(k).rank, refined_front(k).tof_total_days, refined_front(k).dv_total_ms);
end





