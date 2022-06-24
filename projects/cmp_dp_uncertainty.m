% This script computes:

% 1. Pick MAP parameters for HMM and SC from Taiyi Wang's 2021 
% JGR paper, run time dependent inversion for pressure changes between 
% 08/15/2018 and 12/22/2020.

% 2. compute distribution of dynamic pressure changes in Halemaumau
% and South Caldera reservoirs between 08/15/2018 and 12/22/2020, 
% accounting for uncertainties in HMM volume.

% Taiyi Wang 2020-2021

%% load distributions
current_dir  = pwd;
idcs   = strfind(current_dir ,'/');
above_dir = current_dir(1:idcs(end)-1); % one level above current directory

addpath(fullfile(above_dir,'/projects/'));
addpath(fullfile(above_dir,'/source_code/'));
addpath(fullfile(above_dir,'/intermediate/gps_loc_and_offset'));
addpath(fullfile(above_dir,'/intermediate/PDF_from_static_inversion'))

load('static_PDF.mat')
load('HMM_V_distr_kyle.mat')
load('gps_offsets.mat', 'offset_20191201_20201215');
load('gps_locs.mat')

% Set up fixed parameter values
mu = 3e9; nu =0.25; V_HMM = 3.94e9; V_SC = 2.5e9;

% reorganize gps data into a vector
data = offset_20191201_20201215;
d = [data(:,1);data(:,2);data(:,3)];

% load displacement time series-------
dike_corr_flag = 1; % whether to remove the dike intrusion displacement signal on 12/01/2020
CALS_flag = 0; % 1 for including CALS: note that the discontinuation of CALS after 11/22/2020 causes a jump in inferred pressure

% load time series data
addpath(fullfile(above_dir,'/intermediate/gps_time_series'))
load('time_series.mat')

if CALS_flag == 0
    disp_series_ENU(:, 2*3+1:3*3) = []; % take out CALS
    gps_sites(3) = [];
    gps_llh(:, 3) = [];
    gps_xy(:, 3) = [];
end

Ndates = length(time_series_dates);
Ncomps = length(disp_series_ENU(1,:));
Nsites = Ncomps/3;

% percentage interval to extract:
CI = 0.95;
CI_l = (1-CI)/2;
CI_U = 1-CI_l;

% reorganize displacement into an array of [E, E, E, ..., N, N, N, ..., U, U, U]
disp_series = zeros(Ndates, Ncomps);
disp_series_temporary = disp_series_ENU;
disp_series(1:end, 1:Nsites) = disp_series_temporary(:,1:3:Ncomps-2); % East
disp_series(1:end, Nsites+1:2*Nsites) = disp_series_temporary(:,2:3:Ncomps-1); % North
disp_series(1:end, 2*Nsites+1:3*Nsites) = disp_series_temporary(:,3:3:Ncomps); % Up

% Compute displacement relative to day 1
disp_day1 = disp_series(1,:);
for i = 1:Ndates
    disp_series(i, :) = disp_series(i, :) - disp_day1;
end

%PS Correct for Dec 2 2020 dike intrusion. Take difference of one day post
%intrusion to one day pre intrusion.
if dike_corr_flag == 1
    idx1 = find(time_series_dates == datetime('2020-12-02 02:00:00'));
    u_dike = disp_series(idx1+2,:) - disp_series(idx1,:);
    disp_series_corr = disp_series;  
    disp_series_corr(idx1+2:end,:) = disp_series_corr(idx1+2:end,:)- u_dike;
    disp_series_corr(idx1+1,:) = disp_series_corr(idx1+2,:);
    
    % use the corrected time series for inversion
    disp_series = disp_series_corr;
end


%% 1. Pick one set of reservoir geometry parameters, fix volumes, and run inversion
HMM_distr_flag = 0; % whether to sample Kyle's HMM volume distribution; if not, fix at the median value

% use MAP values from Wang et al., 2021
param_MAP = [460, 350, 2180, 1.78, 1890, -3030, 3630, 0.14, 117, -44];

params = param_MAP;

param_labels = ["dx HMM (m)", "dy HMM (m)", "dz HMM (m)", "aspect ratio HMM", "dx SC (m)", "dy SC (m)", "dz SC(m)", "aspect ratio SC", "dip SC", "strike SC"];

% Invert for the best fit pressure for each time step
figure;
for i = 1:10
    param_distr = static_inversion_PDF(:, i);
    subplot(2, 5, i)
    histogram(param_distr, 100, 'EdgeColor', 'none');
    xline(params(i), 'r', 'LineWidth', 1);
    xlabel(param_labels(i));
end

HMM_series = zeros(Ndates, 1);
SC_series = zeros(Ndates, 1);

for i = 1:Ndates
    disp = disp_series(i,:)'; % measured displacement on this day relative to 
    
    min_V = min(V);
    postCol_V = V'- 0.8e9;
    postCol_V = postCol_V(postCol_V > min_V);
        
    GF = cmp_GFs([params, NaN], gps_xy, mu, nu, V_HMM, V_SC, HMM_distr_flag);
    [du, m] = cmp_lsq(GF{1, 1}, disp);
       
    % Store pressure change
    HMM_series(i) = m(1);
    SC_series(i) = m(2);
end

figure;
plot(time_series_dates, HMM_series./1e6, 'r-', 'LineWidth', 2); hold on;
plot(time_series_dates, SC_series./1e6, 'b-', 'LineWidth', 2);
xlabel('days'); ylabel('pressure (MPa)');
title('V_{HMM} = 3.94 km^3, V_{SC} = 2.5 km^3, other params at MAP')
legend('HMM', 'SC')


%% 2. Compute the distribution of dynamic pressure changes for 08/15/18 - 12/22/20
HMM_distr_flag = 1; % whether to sample Kyle's HMM volume distribution; if not, fix at the median value

% Compute static displacement, organized by Nstation X Ncomp (08/15/2018 - 11/22/2020)
du_data = zeros(Nsites, 3);

for i = 1:Nsites
    dE = disp_series(end, i);
    dN = disp_series(end, i+Nsites);
    dZ = disp_series(end, i+2*Nsites);
    
    du_data(i, :) = [dE, dN, dZ];
end

% Invert for the best fit pressure for each time step
Nsamples = 1e3;
m_distr_series = cell(Ndates, 1);
du_distr_series = cell(Ndates, 1);
HMM_CI_series = zeros(Ndates, 2);
SC_CI_series = zeros(Ndates, 2);
HMM_median_series = zeros(Ndates, 1);
SC_median_series = zeros(Ndates, 1);
GF_distr_series = cell(Ndates, 1);

for i = 1:Ndates
    disp = disp_series(i,:)'; % measured displacement on this day relative to beginning of time series

    samples = sample_PDF(static_inversion_PDF, Nsamples, 'true');
        
    min_V = min(V);
    postCol_V = V'- 0.8e9;
    postCol_V = postCol_V(postCol_V > min_V);
        
    HMM_V_sample = sample_PDF(postCol_V, Nsamples, 'true'); % sample Kyle's HMM volume distributions minus co-collapse volume reduction
    GFs = cmp_GFs([samples, HMM_V_sample], gps_xy, mu, nu, V_HMM, V_SC, HMM_distr_flag);
    GF_distr_series{i} = GFs;
        
    % preallocation 
    m_distr = zeros(2, Nsamples);       % pressure change distributions [HMM; SC]
    du_distr = zeros(Ncomps, Nsamples); % displacement distributions [E; N; U]

    for j = 1:Nsamples
        [du, m] = cmp_lsq(GFs{j}, disp);
        m_distr(:,j) = real(m); % the GFs computed from spheroid occasionally has a very small, or zero complex component
        du_distr(:,j) = du;
    end
    
    % Store distributions of pressure and displacement
    m_distr_series{i} = m_distr;
    du_distr_series{i} = du_distr;
    
    % compute confidence interval
    [f1,xi1] = ksdensity(m_distr(1,:)); % HMM
    [f2,xi2] = ksdensity(m_distr(2,:)); % SC
    
    cumProb1 = cumtrapz(xi1, f1); % cumulative probability
    [~,idx_l1] = min(abs(cumProb1 - CI_l)); 
    [~,idx_m1] = min(abs(cumProb1 - 0.5)); 
    [~,idx_U1] = min(abs(cumProb1 - CI_U)); 
    l1 = xi1(idx_l1);
    U1 = xi1(idx_U1);
    median1= xi1(idx_m1);
        
    cumProb2 = cumtrapz(xi2, f2); % cumulative probability
    [~,idx_l2] = min(abs(cumProb2 - CI_l)); 
    [~,idx_m2] = min(abs(cumProb2 - 0.5));
    [~,idx_U2] = min(abs(cumProb2 - CI_U)); 
    l2 = xi2(idx_l2);
    U2 = xi2(idx_U2);
    median2= xi2(idx_m2);
        
    
    % Store confidence interval for this time step
    HMM_CI_series(i, :) = [l1, U1];
    SC_CI_series(i, :) = [l2, U2];
        
    % Store median
    % for pressure change
    HMM_median_series(i) = median1;
    SC_median_series(i) = median2;      
end


%% Plot time series inversion results

% 1. Plot comparison with measured static displacement, with surface
% projection of spheroids

% find the predicted displacement corresponding to median net pressure
% estimates for HMM in the post-collapse period
dp_end = real(m_distr_series{end}); 

% compute sample index of pressure medians
[f1,xi1] = ksdensity(dp_end(1,:)); % HMM
    
cumProb1 = cumtrapz(xi1, f1); % cumulative probability
[~,idx_m1] = min(abs(cumProb1 - 0.5)); 

du_pred_end = du_distr_series{end};
du_pred = du_pred_end(:, idx_m1); % displacement corresponding to median pressure change 
du_pred = [du_pred(1:Nsites), du_pred(Nsites+1:Nsites*2), du_pred(2*Nsites+1:Nsites*3)];


% Sort measured and predicted vertical components into positive and
% negative
data_pos_idx = find(du_data(:,3)>=0);
data_neg_idx = find(du_data(:,3)<0);
pred_pos_idx = find(du_pred(:,3)>=0);
pred_neg_idx = find(du_pred(:,3)<0);

% obtain the geometry parameters for HMM and SC associated with the median
% HMM inverted pressure values at the last time step
s = samples(idx_m1, :); % geometry parameters associated with median HMM pressure change

dx_HMM = s(1); dy_HMM = s(2); dz_HMM = -s(3); alpha_HMM = s(4); 
dx_SC = s(5); dy_SC = s(6); dz_SC = -s(7); alpha_SC = s(8); 
dip_SC = s(9); strike_SC = s(10); 

V_HMM = HMM_V_sample(idx_m1);

a_HMM = (3*V_HMM/(4*pi)*alpha_HMM^2)^(1/3);
b_HMM = a_HMM/alpha_HMM;

a_SC = (3*V_SC/(4*pi)*alpha_SC^2)^(1/3);
b_SC = a_SC/alpha_SC;

% input for spheroid code
m_HMM = [a_HMM; b_HMM; 90; 0; dx_HMM; dy_HMM; dz_HMM; 1e6];
m_SC = [a_SC; b_SC; dip_SC; strike_SC; dx_SC; dy_SC; dz_SC; 1e6];

load('Hawaii_topo.mat')

cScale = 5;
qScale = 5;

figure;
plot(line_xy(1,:), line_xy(2,:), 'k-'); hold on;
quiver(gps_xy(1,:)./1000, gps_xy(2,:)./1000, du_data(:,1)'.*qScale, du_data(:,2)'.*qScale, 'k', 'LineWidth', 2, 'AutoScale', 'off');
quiver(gps_xy(1,:)./1000, gps_xy(2,:)./1000, du_pred(:,1)'.*qScale, du_pred(:,2)'.*qScale, 'r', 'LineWidth', 2, 'AutoScale', 'off');
quiver(6, 4, -0.1.*qScale, 0.*qScale, 'k', 'LineWidth', 2, 'AutoScale', 'off');
viscircles([gps_xy(1,data_pos_idx)'./1000, gps_xy(2,data_pos_idx)'./1000], abs(du_data(data_pos_idx, 3))*cScale, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2);
viscircles([gps_xy(1,data_neg_idx)'./1000, gps_xy(2,data_neg_idx)'./1000], abs(du_data(data_neg_idx, 3))*cScale, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
viscircles([gps_xy(1,pred_pos_idx)'./1000, gps_xy(2,pred_pos_idx)'./1000], abs(du_pred(pred_pos_idx, 3))*cScale , 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2);
viscircles([gps_xy(1,pred_neg_idx)'./1000, gps_xy(2,pred_neg_idx)'./1000], abs(du_pred(pred_neg_idx, 3))*cScale , 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
viscircles([6, 4], 0.1*cScale, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 2); % legend

mdl = [m_HMM,m_SC];

for i = 1:2
    m = mdl(:,i);
    a = m(1);
    b = m(2);
    dip = mod(m(3),180);
    ct = cosd(dip);
    st = sind(dip);
    strike = 90 - m(4);
    cs = cosd(strike);
    ss = sind(strike);
    X1 = m(5);
    X2 = m(6);
    X3 = m(7);
    strength = m(8);
    
    FONTNAME = 'Consolas';
    FONTSIZE = 12;
    R1 = [st 0 ct;0 1 0; ct 0 -st];
    R2 = [cs -ss 0;ss cs 0; 0  0 1];

    n = 50;
    [X,Y,Z] = ellipsoid(0,0,0,b,b,a,n);

    coords = [X(:) Y(:) Z(:)]';
    epts = bsxfun(@plus,R2*R1*coords,[X1;X2;X3]);
    nhat = R2*R1*bsxfun(@rdivide,2*coords,[b^2;b^2;a^2]);
    nhat = bsxfun(@rdivide,nhat,sqrt(sum(nhat.^2)));

    X(:) = epts(1,:)/1000;
    Y(:) = epts(2,:)/1000;
    Z(:) = epts(3,:)/1000; % relative to surface
    
    [sz1, sz2] = size(X);
    RGB(:,:,1) = repelem(0.5, sz1, sz2);
    RGB(:,:,2) = repelem(0.5, sz1, sz2);
    RGB(:,:,3) = repelem(0.5, sz1, sz2);
    hs = surf(X,Y,Z, RGB, 'EdgeColor', 'none');
    
    axis equal
    set(gca,'FontSize',FONTSIZE,'FontName',FONTNAME)
    L = a/1000+1;
    axis([X1/1000+[-L L] X2/1000+[-L L] X3/1000+[-L L]])

    shading interp
    lightangle(40,-15)
    hs.FaceLighting = 'gouraud';
    hs.AmbientStrength = 0.3;
    hs.DiffuseStrength = 0.8;
    hs.SpecularStrength = 0.9;
    hs.SpecularExponent = 25;
    hs.BackFaceLighting = 'unlit';
    
    view([0, 90]);
    grid on;
    hold on;
end

text(5.5, 3.5, '10 cm', 'FontSize', 15);
text(gps_xy(1,:)./1000, gps_xy(2,:)./1000, gps_sites, 'FontSize', 15)
xlim([-6, 8]);ylim([-8, 6]);
xlabel('East (km)', 'FontSize', 15);ylabel('North (km)', 'FontSize', 15)
title('predicted displacement from median HMM dp, 20180815-20201222')
set(gca,'FontSize', 15)




