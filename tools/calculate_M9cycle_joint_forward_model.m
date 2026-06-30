function model = calculate_M9cycle_joint_forward_model(radii, ring_taus, cfg)
%CALCULATE_M9CYCLE_JOINT_FORWARD_MODEL Evaluate predictions for one MCMC state.
%
% This function contains the mechanical forward calculation for the 1998 and
% 2009 velocity fields. The main MCMC script handles sampling, constraints,
% bookkeeping, and file output.

radii_1 = radii(1:end/2);
radii_2 = radii(1+end/2:end);
ring_taus_1 = ring_taus(1:end/2);
ring_taus_2 = ring_taus(1+end/2:end);

Ravg = average_asperity_radii(cfg.asp_points);

i_locked_1 = get_locked_indices_asperities( ...
    radii_1, cfg.asp_centers, cfg.asp_points, cfg.asp_xy, ...
    cfg.Rupture_Scale, cfg.centerx, cfg.centery);
ring_tau_1 = get_ringtau_asperities( ...
    ring_taus_1, cfg.D./Ravg, radii, cfg.asp_centers, cfg.asp_points, ...
    cfg.asp_xy, cfg.Rupture_Scale, cfg.centerx, cfg.centery, i_locked_1);

i_locked_2 = get_locked_indices_asperities( ...
    radii_2, cfg.asp_centers, cfg.asp_points, cfg.asp_xy, ...
    cfg.Rupture_Scale, cfg.centerx, cfg.centery);
ring_tau_2 = get_ringtau_asperities( ...
    ring_taus_2, cfg.D./Ravg, radii, cfg.asp_centers, cfg.asp_points, ...
    cfg.asp_xy, cfg.Rupture_Scale, cfg.centerx, cfg.centery, i_locked_2);

[creep_rate_1, gmres_flag_1, inc_slip_1] = solve_creep_rate(i_locked_1, ring_tau_1, cfg);
[creep_rate_2, gmres_flag_2] = solve_creep_rate(i_locked_2, ring_tau_2, cfg, inc_slip_1);

% Estimate the scale of the viscoelastic contribution for this proposal.
ddhat_no_visc = [cfg.GG*(cfg.rates-creep_rate_1); ...
                 cfg.GG*(cfg.rates-creep_rate_2)] + [cfg.Vpost1; cfg.Vpost2]/3;
residual_no_visc = cfg.dd - ddhat_no_visc;
visco_scale = [cfg.Vvisc/3; cfg.Vvisc/3]\residual_no_visc;

dhat_1 = cfg.G*(cfg.rates-creep_rate_1) + cfg.Vpost1 + visco_scale*cfg.Vvisc;
ddhat_1 = cfg.GG*(cfg.rates-creep_rate_1) + cfg.Vpost1/3 + visco_scale*cfg.Vvisc/3;

dhat_2 = cfg.G*(cfg.rates-creep_rate_2) + cfg.Vpost2 + visco_scale*cfg.Vvisc;
ddhat_2 = cfg.GG*(cfg.rates-creep_rate_2) + cfg.Vpost2/3 + visco_scale*cfg.Vvisc/3;

dhat = [dhat_1; dhat_2];
ddhat = [ddhat_1; ddhat_2];
residual = cfg.dd - ddhat;

model.logrho = -0.5*(residual'*residual);
model.dhat = dhat;
model.dhat_1 = dhat_1;
model.dhat_2 = dhat_2;
model.ddhat = ddhat;
model.i_locked = [i_locked_1; i_locked_2];
model.i_locked_1 = i_locked_1;
model.i_locked_2 = i_locked_2;
model.ring_tau_1 = ring_tau_1;
model.ring_tau_2 = ring_tau_2;
model.creep_rate_1 = creep_rate_1;
model.creep_rate_2 = creep_rate_2;
model.visco_scale = visco_scale;
model.gmres_flag_1 = gmres_flag_1;
model.gmres_flag_2 = gmres_flag_2;

end

function Ravg = average_asperity_radii(asp_points)
Ravg = zeros(1, length(asp_points));

for k = 1:length(asp_points)
    Ravg(k) = mean(sqrt(asp_points{k}(:,1).^2 + asp_points{k}(:,2).^2));
end

end

function [creep_rate, flag, inc_slip] = solve_creep_rate(asp_index, ring_tau, cfg, x0)
bslip = zeros(size(cfg.el, 1), 1);
bslip(asp_index) = cfg.srate(asp_index);

dtau_rate = cfg.scale*hm_mvp('mvp', cfg.gamb.hmat.id, -bslip);
dtau = -dtau_rate*cfg.dT + ring_tau;

rs = cfg.num(~asp_index);
rhs = dtau(~asp_index)/cfg.scale;

if nargin < 4
    [s, flag] = gmres( ...
        @(x)mvp(x, cfg.gamb.hmat.id, rs, rs, ~asp_index), ...
        rhs, [], cfg.tol, cfg.maxit);
else
    [s, flag] = gmres( ...
        @(x)mvp(x, cfg.gamb.hmat.id, rs, rs, ~asp_index), ...
        rhs, [], cfg.tol, cfg.maxit, [], [], x0(~asp_index));
end

inc_slip = zeros(size(cfg.el, 1), 1);
inc_slip(~asp_index) = s;

U = cfg.srate*cfg.dT + inc_slip;
U(asp_index) = 0;

creep_rate = U/cfg.dT*1000;

end
