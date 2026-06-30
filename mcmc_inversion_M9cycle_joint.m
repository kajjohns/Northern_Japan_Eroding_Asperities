%MCMC_INVERSION_M9CYCLE_JOINT Run the Northern Japan M9-cycle MCMC inversion.
%
% Run setup_NorthernJapan.m first. This script assumes the setup variables are
% already in the workspace, then samples asperity radius scales and process-zone
% ring stresses for the 1998 and 2009 velocity fields.

%% User Configuration
% Scale the computed Sanriku postseismic velocities. The original model
% amplitudes are reduced before being added to the elastic predictions.
Ve_post1 = 0.75*Ve_Sanriku_1998;
Vn_post1 = 0.75*Vn_Sanriku_1998;

Ve_post2 = 0.75*Ve_Sanriku_2009;
Vn_post2 = 0.75*Vn_Sanriku_2009;

% Starting asperity radius scale for 1998 and 2009.
radii = [0.7*ones(length(asp_points),1); 0.5*ones(length(asp_points),1)];
stepsize_radii = 0.4*ones(size(radii));

% Starting accumulated stress on the process-zone ring around asperities (Pa).
ring_taus = [0.15e6*ones(length(asp_points),1); 0.15e6*ones(length(asp_points),1)];
stepsize_ring_tau = 1.5e5*ones(size(radii));

% continuing=true resumes from folder_name. continuing=false overwrites the
% output text files in folder_name.
continuing = true;

D = 10;      % Process-zone ring width, km.
mu = 3e10;   % Elastic shear modulus, Pa.

% Optional alternate folder for initial parameter values.
startval = false;
starting_folder = 'NoJapan_outputs_M9cycle_joint';

% Folder for MCMC output files.
folder_name = 'NoJapan_outputs_M9cycle_joint';

%% Static Model Setup
% Fully relaxed M9-cycle response used as a viscoelastic basis function.
load(fullfile('data', 'visco_cycle_M9_contribution.mat'))

addpath ./hmmvp0.16
addpath ./tools

[gamb.hmat.id, nnz] = hm_mvp('init', gamb.hmat.savefn, 4);

patch_stuff = make_triangular_patch_stuff(el, nd);
centroids = patch_stuff.centroids_faces;
strikes = patch_stuff.strike_faces;

% Convert plate rate from m/yr to mm/yr for comparison to GPS velocities.
srate = vp;
rates = srate*1000;

% Horizontal Green's function matrices. GG applies the assumed 3 mm/yr data
% uncertainty used in the likelihood calculation.
G = [Ge; Gn];
GG = [Ge/3; Gn/3];
dd = [Ve1/3; Vn1/3; Ve2/3; Vn2/3];

num = 1:size(el,1);
maxit = 200;
tol = 1e-4;

% Stress Green's functions were computed with km and mu=1; scale to Pa.
scale = mu*10^-3;
dT = 10;  % Time step for process-zone propagation, years.

%% Resume or Seed the Chain
if continuing
    load(fullfile(folder_name, 'M_radii.txt'))
    radii = M_radii(end,:)';

    load(fullfile(folder_name, 'M_ring_tau.txt'))
    ring_taus = M_ring_tau(end,:)';
end

if startval
    load(fullfile(starting_folder, 'M_radii.txt'))
    radii = M_radii(end,:)';

    load(fullfile(starting_folder, 'M_ring_tau.txt'))
    ring_taus = M_ring_tau(end,:)';
end

X = [radii; ring_taus];
stepsize = [stepsize_radii; stepsize_ring_tau];

Vvisc = [Ve_inter_visc; Vn_inter_visc];
Vpost1 = [Ve_post1; Vn_post1];
Vpost2 = [Ve_post2; Vn_post2];

%% Forward-Model Configuration
% Keep the forward model inputs in one struct so each MCMC proposal evaluates
% the same physical configuration with only radii and ring_taus changing.
forward_cfg = struct();
forward_cfg.asp_centers = asp_centers;
forward_cfg.asp_points = asp_points;
forward_cfg.asp_xy = asp_xy;
forward_cfg.Rupture_Scale = Rupture_Scale;
forward_cfg.centerx = centerx;
forward_cfg.centery = centery;
forward_cfg.D = D;
forward_cfg.el = el;
forward_cfg.srate = srate;
forward_cfg.rates = rates;
forward_cfg.scale = scale;
forward_cfg.gamb = gamb;
forward_cfg.num = num;
forward_cfg.tol = tol;
forward_cfg.maxit = maxit;
forward_cfg.dT = dT;
forward_cfg.G = G;
forward_cfg.GG = GG;
forward_cfg.dd = dd;
forward_cfg.Vpost1 = Vpost1;
forward_cfg.Vpost2 = Vpost2;
forward_cfg.Vvisc = Vvisc;

forward_model = calculate_M9cycle_joint_forward_model(radii, ring_taus, forward_cfg);
logrho = forward_model.logrho;
dhat = forward_model.dhat;
i_locked = forward_model.i_locked;
creep_rate_1 = forward_model.creep_rate_1;
creep_rate_2 = forward_model.creep_rate_2;
visco_scale = forward_model.visco_scale;

%% Chain State and Output Files
Xprev = X;
logrhoprev = logrho;
dhatprev = dhat;
i_locked_prev = i_locked;
radii_prev = radii;
ring_tau_prev = ring_taus;
creep_rate_prev_1 = creep_rate_1;
creep_rate_prev_2 = creep_rate_2;
visco_scale_prev = visco_scale;

if ~continuing
    fid = fopen(fullfile(folder_name, 'M_radii.txt'), 'w'); fclose(fid);
    fid = fopen(fullfile(folder_name, 'M_ring_tau.txt'), 'w'); fclose(fid);
    fid = fopen(fullfile(folder_name, 'logrho.txt'), 'w'); fclose(fid);
    fid = fopen(fullfile(folder_name, 'dhat.txt'), 'w'); fclose(fid);
    fid = fopen(fullfile(folder_name, 'locked_index.txt'), 'w'); fclose(fid);
    fid = fopen(fullfile(folder_name, 'creep_rates.txt'), 'w'); fclose(fid);
    fid = fopen(fullfile(folder_name, 'visco_scale.txt'), 'w'); fclose(fid);
end

fidM_radii = fopen(fullfile(folder_name, 'M_radii.txt'), 'a');
fidM_tau = fopen(fullfile(folder_name, 'M_ring_tau.txt'), 'a');
fidrho = fopen(fullfile(folder_name, 'logrho.txt'), 'a');
fiddhat = fopen(fullfile(folder_name, 'dhat.txt'), 'a');
fid_locked = fopen(fullfile(folder_name, 'locked_index.txt'), 'a');
fid_creeping = fopen(fullfile(folder_name, 'creep_rates.txt'), 'a');
fid_scale = fopen(fullfile(folder_name, 'visco_scale.txt'), 'a');

%% MCMC Sampling Loop
numparams = length(X);
numaccept = 0;

rand('state', sum(100*clock));
for iter = 1:10^6

    % Step through parameters one at a time, perturbing a single element.
    count = mod(iter, numparams) + 1;
    r = (-1)^round(rand(1))*rand(1);
    r = r*stepsize(count);
    X(count) = X(count) + r;

    radii = X(1:length(radii));
    ring_taus = X(1+length(radii):end);

    radii_1 = radii(1:end/2);
    radii_2 = radii(1+end/2:end);
    ring_taus_1 = ring_taus(1:end/2);
    ring_taus_2 = ring_taus(1+end/2:end);

    % Reject unphysical proposals before the expensive forward calculation.
    valid_proposal = sum(radii_1 < 0.5 | radii_1 < radii_2 | ...
        ring_taus_1 > ring_taus_2 | ring_taus_1 < 0 | ring_taus_2 < 0 | ...
        ring_taus_1 > 7e5 | ring_taus_2 > 7e5 | ...
        radii_1 < 0 | radii_2 < 0 | radii_1 > 1 | radii_2 > 1) == 0;

    if valid_proposal
        proposal = calculate_M9cycle_joint_forward_model(radii, ring_taus, forward_cfg);
        logrho2 = proposal.logrho;
        accept = metropolis_log(1, logrho, 1, logrho2);
    else
        accept = 0;
    end

    if accept == 1
        logrho = logrho2;

        Xprev = X;
        dhat = proposal.dhat;
        dhatprev = dhat;
        logrhoprev = logrho;

        i_locked = proposal.i_locked;
        i_locked_prev = i_locked;
        radii_prev = radii;
        ring_tau_prev = ring_taus;

        creep_rate_1 = proposal.creep_rate_1;
        creep_rate_2 = proposal.creep_rate_2;
        creep_rate_prev_1 = creep_rate_1;
        creep_rate_prev_2 = creep_rate_2;

        visco_scale = proposal.visco_scale;
        visco_scale_prev = visco_scale;

        numaccept = numaccept + 1;
    else
        X = Xprev;
        dhat = dhatprev;
        logrho = logrhoprev;
        i_locked = i_locked_prev;
        radii = radii_prev;
        ring_taus = ring_tau_prev;
        creep_rate_1 = creep_rate_prev_1;
        creep_rate_2 = creep_rate_prev_2;
        visco_scale = visco_scale_prev;
    end

    % Write one row after every full sweep through the parameter vector.
    if count == 1
        fprintf(fidM_radii, '\n', ' '); fprintf(fidM_radii, '%6.8f\t', radii');
        fprintf(fidrho, '\n', ' '); fprintf(fidrho, '%6.8f\t', logrho);
        fprintf(fiddhat, '\n', ' '); fprintf(fiddhat, '%6.8f\t', dhat');
        fprintf(fidM_tau, '\n', ' '); fprintf(fidM_tau, '%6.8f\t', ring_taus');
        fprintf(fid_locked, '\n', ' '); fprintf(fid_locked, '%6.8f\t', i_locked');
        fprintf(fid_creeping, '\n', ' '); fprintf(fid_creeping, '%6.8f\t', [creep_rate_1; creep_rate_2]');
        fprintf(fid_scale, '\n', ' '); fprintf(fid_scale, '%6.8f\t', visco_scale);

        disp(['Completed sample number ' num2str(iter)])
        disp(['Acceptance rate: ' num2str(numaccept/iter*100)])
        disp(['Visco scale: ' num2str(visco_scale)])
        disp(['Log likelihood: ' num2str(logrho)])
    end

    Logrhos(iter) = logrho;

    figure(1)
    plot(Logrhos)
    drawnow
end
