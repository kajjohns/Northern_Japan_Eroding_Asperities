%PLOT_FANCY_MCMC_INVERSION_JOINT Make summary figures from MCMC output.
%
% Run setup_NorthernJapan.m first so the mesh, asperity geometry, velocity
% fields, and plotting helpers are available in the workspace. This script
% then loads the text output from mcmc_inversion_M9cycle_joint.m and makes
% maps and diagnostic plots for coupling, locking, process-zone stressing,
% velocity fits, moment rates, and asperity shrinkage.

addpath ./tools
addpath ./tools/dem
addpath ./tools/ne_10m_coastline

if ~exist('asp_points', 'var') || ~exist('xysites', 'var') || ~exist('D', 'var')
    error('Run setup_NorthernJapan.m before plot_fancy_MCMC_inversion_joint.m.')
end

%% User Options
folder_name = 'NoJapan_outputs_M9cycle_joint';
discard = 0;  % Number of initial MCMC samples to discard as burn-in.

%% Load MCMC Output
load(fullfile(folder_name, 'logrho.txt'))
figure
plot(logrho)
title('MCMC Log Likelihood')

disp('Loading result files...')

load(fullfile(folder_name, 'M_radii.txt'))
load(fullfile(folder_name, 'M_ring_tau.txt'))
load(fullfile(folder_name, 'dhat.txt'))
load(fullfile(folder_name, 'locked_index.txt'))
load(fullfile(folder_name, 'creep_rates.txt'))
load(fullfile(folder_name, 'visco_scale.txt'))

% Remove burn-in samples before computing posterior summaries.
M_radii(1:discard,:) = [];
M_ring_tau(1:discard,:) = [];
dhat(1:discard,:) = [];
locked_index(1:discard,:) = [];
creep_rates(1:discard,:) = [];
visco_scale(1:discard) = [];

M_radii_1 = M_radii(:,1:end/2);
M_radii_2 = M_radii(:,1+end/2:end);

M_ring_tau_1 = M_ring_tau(:,1:end/2);
M_ring_tau_2 = M_ring_tau(:,1+end/2:end);


locked_index_1 = locked_index(:,1:end/2);
locked_index_2 = locked_index(:,1+end/2:end);

creep_rates_1 = creep_rates(:,1:end/2);
creep_rates_2 = creep_rates(:,1+end/2:end);

dhat_1 = dhat(:,1:end/2);
dhat_2 = dhat(:,1+end/2:end);

%% Reconstruct Process-Zone Stressing Rates
disp('Reconstructing stressing rates...')

clear Ring_Taus_1
clear Ring_Taus_2

% Ravg normalizes the process-zone width for each polygon asperity. This
% matches the calculation inside calculate_M9cycle_joint_forward_model.m.
Ravg = zeros(1, length(asp_points));
for k = 1:length(asp_points)
    Ravg(k) = mean(sqrt(asp_points{k}(:,1).^2 + asp_points{k}(:,2).^2));
end

for k=1:size(locked_index_1,1)

    % 1998 posterior sample.
    radii = M_radii_1(k,:)';
    ring_taus = M_ring_tau_1(k,:)';

    i_locked = logical(locked_index_1(k,:)');

    Ring_Taus_1(:,k) = get_ringtau_asperities(ring_taus,D./Ravg,radii,asp_centers,asp_points,asp_xy,Rupture_Scale,centerx,centery,i_locked);


    % 2009 posterior sample.
    radii = M_radii_2(k,:)';
    ring_taus = M_ring_tau_2(k,:)';

    i_locked = logical(locked_index_2(k,:)');

    Ring_Taus_2(:,k) = get_ringtau_asperities(ring_taus,D./Ravg,radii,asp_centers,asp_points,asp_xy,Rupture_Scale,centerx,centery,i_locked);


end

%% Shared Plot Settings
c1 = colormap(flipud(slanCM('inferno')));
c2 = colormap(flipud(slanCM('Blues')));
mycmap = [c2(1:2:end,:); c1(1:2:end,:)];

nd_ll = local2llh(nd(:,1:2)',fliplr(origin))';
sites_llh = local2llh(xysites',fliplr(origin))';

clear alpha
LA = 45;  % Light azimuth for shaded relief maps.

%% Coupling, Locking, and Process-Zone Stressing Maps
lat_range = [35.5 41.5];
lon_range = [138 145];
[A,lats,lons] = plot_global_etopo1(lat_range,lon_range);

bbox = [lon_range(1) lat_range(1); lon_range(2) lat_range(2)];

figure
hfig = tiledlayout(2,3,'TileSpacing','compact');

ax=nexttile(3);
[h,I,z]=dem(lons,lats,A,'LatLon','LandColor',.8*ones(359,3),'SeaColor',.95*ones(359,3),'Zlim',[-2000 2000],'Contrast',1,'Azimuth',LA);

hold on
h=trisurf(el,nd_ll(:,1),nd_ll(:,2),nd(:,3)+100,1-mean(creep_rates_1,1)'./(srate*1000),'edgecolor','none'); 
colorbar
title('1998 Coupling Ratio')
caxis([-1 1])
colormap(ax,mycmap)

col = 'w';
plot_asperities_llh

plot_coast(bbox)
xlim(lon_range)
ylim([lat_range])
alpha(h,0.5)   

ylim(lat_range)
xlim(lon_range)

ax=nexttile(6);
[h,I,z]=dem(lons,lats,A,'LatLon','LandColor',.8*ones(359,3),'SeaColor',.95*ones(359,3),'Zlim',[-2000 2000],'Contrast',1,'Azimuth',LA);

hold on
h=trisurf(el,nd_ll(:,1),nd_ll(:,2),nd(:,3)+100,1-mean(creep_rates_2,1)'./(srate*1000),'edgecolor','none'); 
colorbar
title('2009 Coupling Ratio')
caxis([-1 1])
colormap(ax,mycmap)

col = 'w';
plot_asperities_llh

plot_coast(bbox)
xlim(lon_range)
ylim([lat_range])
alpha(h,0.5)   

ylim(lat_range)
xlim(lon_range)

% Locking probability is the posterior mean of the locked-patch indicator.
ax=nexttile(2);
[h,I,z]=dem(lons,lats,A,'LatLon','LandColor',.8*ones(359,3),'SeaColor',.9*ones(359,3),'Zlim',[-2000 2000],'Contrast',1,'Azimuth',LA);
hold on

h=trisurf(el,nd_ll(:,1),nd_ll(:,2),nd(:,3)+100,mean(locked_index_1),'edgecolor','none'); 
colormap(ax,flipud(slanCM('heat')));
colorbar     
caxis([0 1])
title('1998 Probability of Locking')
plot_coast(bbox)
xlim(lon_range)
ylim([lat_range])
alpha(h,0.5)   

col = 'b';
plot_asperities_llh

ylim(lat_range)
xlim(lon_range)

ax=nexttile(5);
[h,I,z]=dem(lons,lats,A,'LatLon','LandColor',.8*ones(359,3),'SeaColor',.9*ones(359,3),'Zlim',[-2000 2000],'Contrast',1,'Azimuth',LA);
hold on

h=trisurf(el,nd_ll(:,1),nd_ll(:,2),nd(:,3)+100,mean(locked_index_2),'edgecolor','none'); 
colormap(ax,flipud(slanCM('heat')));
colorbar     
caxis([0 1])
title('2009 Probability of Locking')
plot_coast(bbox)
xlim(lon_range)
ylim([lat_range])
alpha(h,0.5)   

col = 'b';
plot_asperities_llh


ylim(lat_range)
xlim(lon_range)
% Process-zone stressing rate reconstructed from posterior ring stresses.
ax=nexttile(1);
[h,I,z]=dem(lons,lats,A,'LatLon','LandColor',.8*ones(359,3),'SeaColor',.9*ones(359,3),'Zlim',[-2000 2000],'Contrast',1,'Azimuth',LA);
hold on

h=trisurf(el,nd_ll(:,1),nd_ll(:,2),nd(:,3)+100,mean(Ring_Taus_1,2)/10,'edgecolor','none'); 
colormap(ax,flipud(slanCM('heat')));
colorbar     
title('1998 Process Zone Stressing Rate')
plot_coast(bbox)
xlim(lon_range)
ylim([lat_range])
alpha(h,0.7)   

ylim(lat_range)
xlim(lon_range)

ax=nexttile(4);
[h,I,z]=dem(lons,lats,A,'LatLon','LandColor',.8*ones(359,3),'SeaColor',.9*ones(359,3),'Zlim',[-2000 2000],'Contrast',1,'Azimuth',LA);
hold on

h=trisurf(el,nd_ll(:,1),nd_ll(:,2),nd(:,3)+100,mean(Ring_Taus_2,2)/10,'edgecolor','none'); 
colormap(ax,flipud(slanCM('heat')));
colorbar     
title('2009 Process Zone Stressing Rate')
plot_coast(bbox)
xlim(lon_range)
ylim([lat_range])
alpha(h,0.7)   

ylim(lat_range)
xlim(lon_range)

set(gcf, 'Position', [31    46   872   546])
set(gcf,'renderer','Painters')

%% Velocity Fit Maps
dhat_e1 = dhat_1(:,1:end/2);
dhat_n1 = dhat_1(:,1+end/2:end);

lat_range = [35.5 41.5];
lon_range = [138 142.5];
[A,lats,lons] = plot_global_etopo1(lat_range,lon_range);

bbox = [lon_range(1) lat_range(1); lon_range(2) lat_range(2)];

figure
hfig2 = tiledlayout(1,3,'TileSpacing','compact');

nexttile(1)
[h,I,z]=dem(lons,lats,A,'LatLon','LandColor',.8*ones(359,3),'SeaColor',.9*ones(359,3),'Zlim',[-2000 2000],'Contrast',1,'Azimuth',LA);
hold on

scale=0.02;
quiver(sites_llh(:,1),sites_llh(:,2),scale*mean(dhat_e1)',scale*mean(dhat_n1)',0,'r','linewidth',1)
quiver(sites_llh(:,1),sites_llh(:,2),scale*Ve1,scale*Vn1,0,'b','linewidth',1)
quiver(sites_llh(:,1),sites_llh(:,2),scale*mean(visco_scale)*Ve_inter_visc,scale*mean(visco_scale)*Vn_inter_visc,0,'k','linewidth',1)

axis equal    
title('1998 Fit to data')
set(gca,'fontsize',15)

quiver(138.5,40,scale*30,0,'k')
text(138.5,40.2,'30 mm/yr')


dhat_e2 = dhat_2(:,1:end/2);
dhat_n2 = dhat_2(:,1+end/2:end);

nexttile(2)
[h,I,z]=dem(lons,lats,A,'LatLon','LandColor',.8*ones(359,3),'SeaColor',.9*ones(359,3),'Zlim',[-2000 2000],'Contrast',1,'Azimuth',LA);
hold on

scale=0.02;
quiver(sites_llh(:,1),sites_llh(:,2),scale*mean(dhat_e2)',scale*mean(dhat_n2)',0,'r','linewidth',1)
quiver(sites_llh(:,1),sites_llh(:,2),scale*Ve2,scale*Vn2,0,'b','linewidth',1)
quiver(sites_llh(:,1),sites_llh(:,2),scale*mean(visco_scale)*Ve_inter_visc,scale*mean(visco_scale)*Vn_inter_visc,0,'k','linewidth',1)
axis equal    
title('2009 Fit to data')
set(gca,'fontsize',15)

quiver(138.5,40,scale*30,0,'k')
text(138.5,40.2,'30 mm/yr')

set(gcf, 'Position', [-38   212   920   374])
set(gcf,'renderer','Painters')

% Simple residual diagnostics using the same 3 mm/yr uncertainty as the
% MCMC likelihood.
resid = [Ve1-mean(dhat_e1)';Vn1-mean(dhat_n1)'];
chi2r_1 = resid'*resid/3^2/length(resid)
var_red_1 = 1-norm(resid)/norm([Ve1;Vn1])

resid = [Ve2-mean(dhat_e2)';Vn2-mean(dhat_n2)'];
chi2r_2 = resid'*resid/3^2/length(resid)
var_red_2 = 1-norm(resid)/norm([Ve2;Vn2])


%% Moment-Rate Summaries
A = patch_stuff.area_faces*10^6;  %convert km^2 to m^2

sr = (srate)-mean(creep_rates_1)'/1000; % m/yr
Mo_total_1998 = 30e9*sum(sr.*A)
Mw_500_total_1998 = (2/3) * log10(500*Mo_total_1998*1e7) - 10.7

locking_prob = mean(locked_index_1)';
Mo_locking_1998 = 30e9*sum(sr.*A.*locking_prob)
Mw_500_locking_1998 = (2/3) * log10(500*Mo_locking_1998*1e7) - 10.7

sr = (srate)-mean(creep_rates_2)'/1000; % m/yr
Mo_total_2009 = 30e9*sum(sr.*A)
Mw_500_total_2009 = (2/3) * log10(500*Mo_total_2009*1e7) - 10.7

locking_prob = mean(locked_index_2)';
Mo_locking_2009 = 30e9*sum(sr.*A.*locking_prob)
Mw_500_locking_2009 = (2/3) * log10(500*Mo_locking_2009*1e7) - 10.7

% Build posterior distributions of total and locked moment rates.
sr = repmat(srate,1,size(creep_rates_1,1))-creep_rates_1'/1000; % m/yr
Mo_total_1998 = 30e9*sum(sr.*repmat(A,1,size(creep_rates_1,1)));

Mo_locking_1998 = 30e9*sum(sr.*repmat(A,1,size(creep_rates_1,1)).*locked_index_1');
Mw_500_locking_1998 = (2/3) * log10(500*Mo_locking_1998*1e7) - 10.7


sr = repmat(srate,1,size(creep_rates_2,1))-creep_rates_2'/1000; % m/yr
Mo_total_2009 = 30e9*sum(sr.*repmat(A,1,size(creep_rates_2,1)));

Mo_locking_2009 = 30e9*sum(sr.*repmat(A,1,size(creep_rates_2,1)).*locked_index_2');
Mw_500_locking_2009 = (2/3) * log10(500*Mo_locking_2009*1e7) - 10.7

figure
histogram(Mo_total_1998,'Normalization','pdf')
hold on
histogram(Mo_total_2009,'Normalization','pdf')
histogram(Mo_locking_1998,'Normalization','pdf')
histogram(Mo_locking_2009,'Normalization','pdf')


xlabel('Moment rate or Mo [N·m/yr]')
ylabel('Counts')
legend('total 1998','total 2009','locked areas 1998','locked areas 2009')
set(gca,'FontSize',15)

% ---- Add second x-axis on top for Mw ----
T = 500;
% Mo is moment *rate* in N·m/yr; convert to 500-yr moment (N·m) then to Mw.
% Using Hanks & Kanamori with dyne-cm: Mw = (2/3)(log10(M0[dyncm]) - 10.7)
Mo2Mw = @(Mo) (2/3) * log10(T * Mo * 1e7) - 10.7;

ax1 = gca;

% Make sure Mo axis only draws ticks at the bottom and no top box line
ax1.XAxisLocation = 'bottom';
ax1.Box = 'off';
ax1.Layer = 'top';    % draw ticks above the bars
drawnow

% Create a transparent axes on top
ax2 = axes('Position', ax1.Position, ...
           'Color', 'none', ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right');

ax2.YLim = ax1.YLim;   % manually sync once
ax2.YTick = [];        % hide y-ticks on the top axis

% Set Mw limits to match Mo limits via the transform
ax2.XLim = Mo2Mw(ax1.XLim);

% Label + styling
xlabel(ax2, 'Equivalent 500-year Mw')
set(ax2,'FontSize',15,'XColor','k')

set(gcf,'renderer','Painters')


%% Evaluate Cattania rule for shrinking (Cattania and Segall, 2019)
% Stress drops are tuned manually for each polygon asperity.
deltaus = 1.5e6*ones(14,1);
deltaus(14)=3.5e6;
deltaus(12)=0.8e6;
deltaus(11)=3.0e6;
deltaus(10)=3.5e6;
deltaus(9)=3.0e6;
deltaus(8)=2.5e6;
deltaus(6)=1e6;
deltaus(5)=1e6;
deltaus(4)=3e6;
deltaus(3)=0.25e6;
deltaus(2)=1.2e6;
deltaus(1)=1.5e6;


mu = 3*10^10;  % Pa.
mup = mu/(1-0.25);

% Rupture radius, converted from km to meters.
R = 1000*sqrt(Rupture_Areas);
cnt=0;
for loop = 1:size(M_radii_1,2)

    cnt=cnt+1;
    if loop==1 || loop==10
        cnt=1;
        figure
        set(gcf, 'Position', [48         295        1071         329])
    end

    subplot(3,3,cnt)
    hold on

    dtau = deltaus(loop);
    mu = 3*10^10;  % Pa.
    
    % Adjust effective shear modulus for geometric factor.
    mup = .5*mu/(1-0.25);


    rc(loop) = pi*mup*srate(1)./(4*dtau);  %m/yr
   

    
    tf = R(loop)./rc(loop);
    ts = linspace(0,max(tf)+10);
   
    rs = sqrt(1-ts.*rc(loop)./R(loop));
    plot(date_mag(loop,1)+ts,real(rs),'linewidth',1)
   
    meanrad1 = mean(M_radii_1(:,loop));
    stdrad1 = std(M_radii_1(:,loop));
    meanrad2 = mean(M_radii_2(:,loop));
    stdrad2 = std(M_radii_2(:,loop));


    errorbar(1998,meanrad1,2*stdrad1,'ro','linewidth',1)
    errorbar(2009,meanrad2,2*stdrad2,'ro','linewidth',1)

    xlabel('calendar year')
    ylabel('radius/R_{rup}')
 
    title(['Asperity' num2str(loop)])
    grid on

    text(date_mag(loop,1),0.2,['\Delta\tau = ' num2str(deltaus(loop)/1e6)])
end



% Posterior shrinkage summaries.
radius_ratio = M_radii_2./M_radii_1;
mean_radius_ratio = mean(radius_ratio)


mean_percent_reduction = mean((M_radii_1-M_radii_2)./M_radii_1);
average_reduction = mean(mean_percent_reduction)

mean_radius_to_rupture_1 = mean(M_radii_1)
mean_radius_to_rupture_2 = mean(M_radii_2)
