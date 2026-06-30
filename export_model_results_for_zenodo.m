%EXPORT_MODEL_RESULTS_FOR_ZENODO Write plotted model results to text files.
%
% Run setup_NorthernJapan.m first so mesh geometry, observation sites, and
% model constants are available in the workspace. This script loads the MCMC
% output files, reconstructs the same posterior quantities shown by
% plot_fancy_MCMC_inversion_joint.m, and writes tab-delimited text files for
% archiving.

addpath ./tools

if ~exist('asp_points', 'var') || ~exist('xysites', 'var') || ~exist('D', 'var')
    error('Run setup_NorthernJapan.m before export_model_results_for_zenodo.m.')
end

%% User Options
folder_name = 'NoJapan_outputs_M9cycle_joint';
output_dir = 'Zenodo_model_results';
discard = 0;  % Number of initial MCMC samples to discard as burn-in.

if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

%% Load MCMC Output
load(fullfile(folder_name, 'M_radii.txt'))
load(fullfile(folder_name, 'M_ring_tau.txt'))
load(fullfile(folder_name, 'dhat.txt'))
load(fullfile(folder_name, 'locked_index.txt'))
load(fullfile(folder_name, 'creep_rates.txt'))

M_radii(1:discard,:) = [];
M_ring_tau(1:discard,:) = [];
dhat(1:discard,:) = [];
locked_index(1:discard,:) = [];
creep_rates(1:discard,:) = [];

M_radii_1 = M_radii(:, 1:end/2);
M_radii_2 = M_radii(:, 1+end/2:end);

M_ring_tau_1 = M_ring_tau(:, 1:end/2);
M_ring_tau_2 = M_ring_tau(:, 1+end/2:end);

locked_index_1 = locked_index(:, 1:end/2);
locked_index_2 = locked_index(:, 1+end/2:end);

creep_rates_1 = creep_rates(:, 1:end/2);
creep_rates_2 = creep_rates(:, 1+end/2:end);

dhat_1 = dhat(:, 1:end/2);
dhat_2 = dhat(:, 1+end/2:end);

%% Mesh Coordinates for Patch-Based Exports
% The plotted trisurf quantities are one value per triangular patch. Each
% output row therefore describes one patch using its centroid plus the lon/lat
% coordinates of its three mesh nodes.
patch_stuff = make_triangular_patch_stuff(el, nd);
centroid_ll = local2llh(patch_stuff.centroids_faces(:, 1:2)', fliplr(origin))';
node_ll = local2llh(nd(:, 1:2)', fliplr(origin))';

%% Reconstruct Process-Zone Stressing Rates
Ravg = zeros(1, length(asp_points));
for k = 1:length(asp_points)
    Ravg(k) = mean(sqrt(asp_points{k}(:,1).^2 + asp_points{k}(:,2).^2));
end

Ring_Taus_1 = zeros(size(el, 1), size(locked_index_1, 1));
Ring_Taus_2 = zeros(size(el, 1), size(locked_index_2, 1));

for k = 1:size(locked_index_1, 1)
    radii = M_radii_1(k,:)';
    ring_taus = M_ring_tau_1(k,:)';
    i_locked = logical(locked_index_1(k,:)');
    Ring_Taus_1(:,k) = get_ringtau_asperities( ...
        ring_taus, D./Ravg, radii, asp_centers, asp_points, asp_xy, ...
        Rupture_Scale, centerx, centery, i_locked);

    radii = M_radii_2(k,:)';
    ring_taus = M_ring_tau_2(k,:)';
    i_locked = logical(locked_index_2(k,:)');
    Ring_Taus_2(:,k) = get_ringtau_asperities( ...
        ring_taus, D./Ravg, radii, asp_centers, asp_points, asp_xy, ...
        Rupture_Scale, centerx, centery, i_locked);
end

%% Write Patch-Based Quantities Plotted with trisurf
coupling_1998 = 1 - bsxfun(@rdivide, creep_rates_1, (srate(:)'*1000));
coupling_2009 = 1 - bsxfun(@rdivide, creep_rates_2, (srate(:)'*1000));

write_patch_file( ...
    fullfile(output_dir, 'coupling_ratio_1998.txt'), ...
    el, centroid_ll, node_ll, mean(coupling_1998, 1)', std(coupling_1998, 0, 1)', ...
    'coupling_ratio');

write_patch_file( ...
    fullfile(output_dir, 'coupling_ratio_2009.txt'), ...
    el, centroid_ll, node_ll, mean(coupling_2009, 1)', std(coupling_2009, 0, 1)', ...
    'coupling_ratio');

write_patch_file( ...
    fullfile(output_dir, 'locking_probability_1998.txt'), ...
    el, centroid_ll, node_ll, mean(locked_index_1, 1)', std(locked_index_1, 0, 1)', ...
    'locking_probability');

write_patch_file( ...
    fullfile(output_dir, 'locking_probability_2009.txt'), ...
    el, centroid_ll, node_ll, mean(locked_index_2, 1)', std(locked_index_2, 0, 1)', ...
    'locking_probability');

write_patch_file( ...
    fullfile(output_dir, 'process_zone_stressing_rate_1998.txt'), ...
    el, centroid_ll, node_ll, mean(Ring_Taus_1, 2)/10, std(Ring_Taus_1, 0, 2)/10, ...
    'process_zone_stressing_rate');

write_patch_file( ...
    fullfile(output_dir, 'process_zone_stressing_rate_2009.txt'), ...
    el, centroid_ll, node_ll, mean(Ring_Taus_2, 2)/10, std(Ring_Taus_2, 0, 2)/10, ...
    'process_zone_stressing_rate');

%% Write Velocity Fit Files
dhat_e1 = dhat_1(:, 1:end/2);
dhat_n1 = dhat_1(:, 1+end/2:end);
dhat_e2 = dhat_2(:, 1:end/2);
dhat_n2 = dhat_2(:, 1+end/2:end);

site_ll = local2llh(xysites', fliplr(origin))';

write_velocity_file( ...
    fullfile(output_dir, 'velocity_fit_1998.txt'), site_ll, ...
    Ve1, Vn1, mean(dhat_e1, 1)', mean(dhat_n1, 1)', ...
    std(dhat_e1, 0, 1)', std(dhat_n1, 0, 1)');

write_velocity_file( ...
    fullfile(output_dir, 'velocity_fit_2009.txt'), site_ll, ...
    Ve2, Vn2, mean(dhat_e2, 1)', mean(dhat_n2, 1)', ...
    std(dhat_e2, 0, 1)', std(dhat_n2, 0, 1)');

disp(['Wrote Zenodo model result files to ' output_dir])

%% Local Helper Functions
function write_patch_file(filename, el, centroid_ll, node_ll, value_mean, value_std, value_name)
node1 = el(:,1);
node2 = el(:,2);
node3 = el(:,3);

T = table( ...
    (1:size(el,1))', ...
    centroid_ll(:,1), centroid_ll(:,2), ...
    node1, node_ll(node1,1), node_ll(node1,2), ...
    node2, node_ll(node2,1), node_ll(node2,2), ...
    node3, node_ll(node3,1), node_ll(node3,2), ...
    value_mean, value_std, ...
    'VariableNames', { ...
        'patch_id', ...
        'centroid_lon', 'centroid_lat', ...
        'node1_id', 'node1_lon', 'node1_lat', ...
        'node2_id', 'node2_lon', 'node2_lat', ...
        'node3_id', 'node3_lon', 'node3_lat', ...
        [value_name '_mean'], [value_name '_std']});

writetable(T, filename, 'Delimiter', '\t', 'FileType', 'text');
end

function write_velocity_file(filename, site_ll, observed_e, observed_n, modeled_e_mean, modeled_n_mean, modeled_e_std, modeled_n_std)
T = table( ...
    (1:size(site_ll,1))', ...
    site_ll(:,1), site_ll(:,2), ...
    observed_e, observed_n, ...
    modeled_e_mean, modeled_n_mean, ...
    modeled_e_std, modeled_n_std, ...
    'VariableNames', { ...
        'site_id', 'lon', 'lat', ...
        'observed_east_mm_per_yr', 'observed_north_mm_per_yr', ...
        'modeled_east_mean_mm_per_yr', 'modeled_north_mean_mm_per_yr', ...
        'modeled_east_std_mm_per_yr', 'modeled_north_std_mm_per_yr'});

writetable(T, filename, 'Delimiter', '\t', 'FileType', 'text');
end
