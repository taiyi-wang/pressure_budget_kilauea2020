% script to determine time to threshold for Kilauea eruption in 2020

current_dir  = pwd;
idcs   = strfind(current_dir ,'/');
above_dir = current_dir(1:idcs(end)-1); % one level above current directory

addpath(fullfile(above_dir,'/output/'))

% load distribution of pressure changes
% Use relative path, assuming in Research/PressureBudget/Surges
load('dp_HMM_12222020.mat');

% convert distribution of HMM pressures to MPa
dp = dp_HMM_12222020*1e-6;   

% load time series of pressure history for nominal case
load('dp_HMM_SC_one_model.mat')

% SOME SMOOTHING VALUES
bw = 50; % bandwidth for ksdensity closer to default value
fltpts = 5;  % points in median filter

% Nominal Pressure Threshold based on known 2020 vent elevation
p_thresh = 11.5;

% Consider a range of elevations and corresponding pressure thresholds
DeltaP_thresh = 5; % in MPa
Nthresh = 11;
pressure_thresh_range = linspace(p_thresh - DeltaP_thresh, ...
    p_thresh + DeltaP_thresh, Nthresh);
DepthRange = 40*(pressure_thresh_range - p_thresh);


%% apply median filter to data
pHmmt = medfilt1(HMM_series*1e-6,fltpts);   % convert to MPa
pSCt = medfilt1(SC_series*1e-6,fltpts);   % convert to MPa

% only consider data to date of eruption in 2020
idx1 = find(time_series_dates == datetime('2020-12-20 02:00:00'));
pHmmt(idx1:end) = [];
pSCt(idx1:end) = [];
time_series_dates(idx1:end) = [];

dy = datenum(time_series_dates);  % convert to decimal year

% make dates relative to Jan1 2018
dy = dy - datenum('01-Jan-2018'); % DoY relative to '01-Jan-2018'

% to test the input
if(1)
    figure; plot(time_series_dates, pHmmt, 'r', time_series_dates, pSCt, 'b', 'LineWidth',3); hold on
    set(gca,'FontSize', 14); hold on
    ylabel('Pressure Change (MPa)', 'FontSize', 16)
    leg = legend('Halemaumau','South Caldera');

    set(leg,'AutoUpdate','off')
    xl = xlim;
    plot(xl, p_thresh*[1 1], 'k--')

end

%% Scale time series corresponding to range of sampled pressure changes
% 
Dpt = zeros(length(pHmmt),length(dp));
for i = 1:length(dp)
    Dpt(:,i) = pHmmt * dp(i)/pHmmt(end);   
end

% test figure
% figure; plot(dy, Dpt)

%% Determine the time at which the threshold is met.
%  For cases where the threshold is never met set the threshold time to a large value.

t_thresh = zeros(length(dp),Nthresh);
for k = 1:Nthresh
    % find index of time when closest to threshold
    Index = zeros(length(dp),1);
    for i = 1:length(dp)
            [~, Index(i)] = min(  abs(Dpt(:,i) - pressure_thresh_range(k)));
    end
    t_thresh(:,k) = dy(Index); % Threshold time in decimal year relative to 2018.0
   
    % Determine cases where the threshold is never met. Two strategies
        %[~, Imax] = max(pHmmt);
        % I = find(Index == Imax);  
        Dptmax = max(Dpt);
        II = find(Dptmax < pressure_thresh_range(k)); % The maximum pressure never reaches threshold
    % It seems both I and II give the same results; maximum pressure criterion is simplest
        t_thresh(II,k) = 2*dy(end); % set threshold time to twice the data length for these cases
end

%% plot CDF of Threshold Time for nominal elevation
% for testing
xi = 390:1:2*max(dy); % dates to consider

k = 6;  % this corresponds to the nominal elevation case
[f,xi,bw] = ksdensity(t_thresh(:,k),xi,'support', 'positive','BoundaryCorrection',...
   'reflection', 'Function','cdf','Bandwidth',bw);

figure; area(xi,f,'facealpha',0.5,'edgecolor','none', 'HandleVisibility','off'); 
hold  on
set(gca,'FontSize', 14); hold on
xlabel('Days From 1/1/2018', 'FontSize', 16)
ylabel('CDF', 'FontSize', 16)
xlim([0 max(dy)])
title('Probability of Reaching Nominal Threshold')


%% Compute CDF for a range of pressure thresholds
xi = 390:1:2*max(dy);
f = zeros(Nthresh,length(xi));

for k = 1:Nthresh
    [f(k,:)] = ksdensity(t_thresh(:,k),xi, 'support', 'positive','BoundaryCorrection',...
    'reflection', 'Function','cdf','Bandwidth',bw);
end

%% Plot for range of thresholds
ti = xi - 2*365+1; % Define time variable relative to Jan 1 2020
[X2,Y2] = meshgrid(ti,DepthRange);

figure;
[MM, cc] = contour(X2,Y2,f,'ShowText','on', 'LineWidth',2);
clabel(MM,cc,'FontSize',14)
set(gca,'FontSize', 14); hold on
grid
xlim([-365 365-12])  % upper limit is day of the eruption
xlabel('DOY 2020', 'FontSize', 16)
ylabel('Depth Relative to 2020 vent, m', 'FontSize', 16)
hold on

ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'),...¨
'YAxisLocation','right','Color','none','YColor','k');
ax2.YLim = [pressure_thresh_range(1) pressure_thresh_range(end)];
ax2.XTick = [];
ylabel(ax2,'Pressure Threshold, MPa', 'FontSize', 16)
set(ax2,'FontSize', 14);
colormap cool


%% interpolate results to get probabilities 100 days before the eruption

Vq = interp2(X2,Y2,f,(365-12-100)*ones(1,length(ti)),Y2(:,1));
disp('Elevation; Probability')
disp([Y2(:,1),Vq(:,1)])
figure; plot(Y2(:,1),Vq(:,1), 'LineWidth',2)

%% Plot CDF for nomina case wrt DOY in 2020

k = 6;
figure; area(ti,f(k,:),'facealpha',0.5,'edgecolor','none', 'HandleVisibility','off'); 
hold  on
set(gca,'FontSize', 14); hold on
xlabel('DOY 2020', 'FontSize', 16)
ylabel('CDF Threshold Exceeded', 'FontSize', 16)

% Day of eruption is 365-12
xlim([-365 365-12])
title(['Pressure Threshold   ',num2str(pressure_thresh_range(k))])

%% inverse cdf to get date at which a threshold probability is met

pi = [0.5, 0.6, 0.7];  %probabilities
k = 6;
t_thresh_crit = ksdensity(t_thresh(:,k),pi,'support', 'positive','BoundaryCorrection',...
    'reflection','Function','icdf','Bandwidth',bw);

% convert to calendar time

crit_date = string(zeros(3));
for i = 1:length(pi)
    crit_date(i,:) = string(doy2cal(t_thresh_crit(i)-2*365+1,2020));    
end

%crit_date_out = [string(pi'), crit_date];
disp('probability,   date')
disp([string(pi'), crit_date])

% as a check final date
%  end_date = doy2cal(dy(end)-2*365+1,2020)

%writematrix(crit_date_out)
