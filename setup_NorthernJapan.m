%SETUP_NORTHERNJAPAN Load shared data and prepare model inputs.
%
% Run this script before mcmc_inversion_M9cycle_joint.m. It loads the mesh,
% Green's functions, observed interseismic velocities, earthquake catalog,
% asperity geometry, and Sanriku postseismic correction fields used by the
% MCMC inversion.

addpath ./hmmvp0.16
addpath ./tools
addpath ./data

%% Mesh, Green's Functions, and Observation Sites
% setup_mesh_minimal.mat provides the triangular mesh, observation sites,
% Green's functions, H-matrix metadata, and origin used throughout the model.
load(fullfile('data', 'setup_mesh_minimal.mat'))
xy = xysites;
dispG = Disp;  % Backslip displacement Green's functions.

% Split interleaved east/north/up displacement rows into horizontal matrices.
Ge = dispG(1:3:end,:);
Gn = dispG(2:3:end,:);

%% Observed Interseismic Velocities
% Columns are east/north velocity for the 1998 and 2009 observation windows.
% Units are mm/yr.
interseismic_velocities = readmatrix(fullfile('data', 'interseismic_velocities_1998_2009.txt'), 'NumHeaderLines', 1);
Ve1 = interseismic_velocities(:, 1);  % 1998 east velocity
Vn1 = interseismic_velocities(:, 2);  % 1998 north velocity
Ve2 = interseismic_velocities(:, 3);  % 2009 east velocity
Vn2 = interseismic_velocities(:, 4);  % 2009 north velocity

% Difference fields are retained for diagnostics and plotting.
DVeast = Ve2 - Ve1;
DVnorth = Vn2 - Vn1;

%% Plate Rate and Earthquake Catalog
npatch = size(el, 1);
vp = 0.085*ones(npatch, 1);  % Plate rate on each patch, m/yr.

earthquake_catalog = readmatrix(fullfile('data', 'earthquake_catalog.txt'), 'NumHeaderLines', 1);
eqdates = earthquake_catalog(:, 1);
eqll = earthquake_catalog(:, 2:3)';
eqmags = earthquake_catalog(:, 4);

eqxy = llh2local([eqll(1,:); eqll(2,:)], [origin(2); origin(1)])';

%% Asperity Geometry and Rupture Scaling
% define_asperities_seismic_M9 returns polygon asperities, circular asperity
% centers, and mesh-centroid coordinates used for locking/ring masks.
[asp_centers, asp_points, asp_xy, centerx, centery, centers] = define_asperities_seismic_M9(origin, centroids);

[date_mag, Rupture_Areas, Rupture_Scale] = make_rupture_dimension_M9(asp_points);

%% Sanriku Postseismic Correction
% Sanriku_sep19_2022.txt is a column export of the original .mat file.
% Acceleration columns are in mm/yr^2; velocity model columns are in m/yr.
sanriku_data = readmatrix(fullfile('data', 'Sanriku_sep19_2022.txt'), 'NumHeaderLines', 1);
xysites_Sanriku = sanriku_data(:, 1:2);
ae_model_Sanriku = sanriku_data(:, 3);
an_model_Sanriku = sanriku_data(:, 4);
au_model_Sanriku = sanriku_data(:, 5);
ve_1998_model_Sanriku = sanriku_data(:, 6);
vn_1998_model_Sanriku = sanriku_data(:, 7);
vu_1998_model_Sanriku = sanriku_data(:, 8);
ve_2009_model_Sanriku = sanriku_data(:, 9);
vn_2009_model_Sanriku = sanriku_data(:, 10);
vu_2009_model_Sanriku = sanriku_data(:, 11);

% Interpolate Sanriku acceleration onto the observation sites. NaNs occur
% outside the interpolation domain and are treated as no correction.
a_e_visc = griddata(xysites_Sanriku(:,1), xysites_Sanriku(:,2), ae_model_Sanriku, xysites(:,1), xysites(:,2));
a_n_visc = griddata(xysites_Sanriku(:,1), xysites_Sanriku(:,2), an_model_Sanriku, xysites(:,1), xysites(:,2));
a_e_visc(isnan(a_e_visc)) = 0;
a_n_visc(isnan(a_n_visc)) = 0;

% Convert acceleration correction to an effective 1998 viscous velocity.
v_e_visc_1998 = (-a_e_visc*11)/1.5;
v_n_visc_1998 = (-a_n_visc*11)/1.5;

% Interpolate modeled Sanriku postseismic velocities onto observation sites
% and convert from m/yr to mm/yr.
Ve_Sanriku_1998 = griddata(xysites_Sanriku(:,1), xysites_Sanriku(:,2), ve_1998_model_Sanriku*1000, xysites(:,1), xysites(:,2), 'nearest');
Vn_Sanriku_1998 = griddata(xysites_Sanriku(:,1), xysites_Sanriku(:,2), vn_1998_model_Sanriku*1000, xysites(:,1), xysites(:,2), 'nearest');

Ve_Sanriku_2009 = griddata(xysites_Sanriku(:,1), xysites_Sanriku(:,2), ve_2009_model_Sanriku*1000, xysites(:,1), xysites(:,2), 'nearest');
Vn_Sanriku_2009 = griddata(xysites_Sanriku(:,1), xysites_Sanriku(:,2), vn_2009_model_Sanriku*1000, xysites(:,1), xysites(:,2), 'nearest');
