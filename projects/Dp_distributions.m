% Script to make probability plots of pressure change @ Kilauea for 2018 to
% December 20 2020. Link to data

current_dir  = pwd;
idcs   = strfind(current_dir ,'/');
above_dir = current_dir(1:idcs(end)-1); % one level above current directory

addpath(fullfile(above_dir,'/source_code/'));
addpath(fullfile(above_dir,'/output/'))
load('dp_HMM_12222020');

% Convert HMM pressure to MPa
dp = dp_HMM_12222020*1e-6;  
med = median(dp); % Median pressure change

% Threshold pressure
pressure_thresh = 11.5;


% if(0)
%     % plot pdf of pressure change as check
%     figure; histogram(dp,'BinMethod','fd','Normalization', 'pdf')
%     set(gca,'FontSize', 14); hold on
%     xlabel('Pressure Change, MPa', 'FontSize', 16)
%     yl = ylim;
%     plot(med*[1 1], yl, 'r','LineWidth',2)
% 
%     % plot cdf of pressure change
%     figure; histogram(dp,'BinMethod','fd','Normalization', 'cdf')
%     set(gca,'FontSize', 14); hold on
%     xlabel('Pressure Change, MPa', 'FontSize', 16)
% end

%% use ksdensity to plot pdf
[f,xi,bw] = ksdensity(dp,'support', 'positive', 'BoundaryCorrection',...
    'reflection');
figure; area(xi,f,'facealpha',0.5,'edgecolor','none', 'HandleVisibility','off'); 
hold  on
set(gca,'FontSize', 14); hold on
xlabel('MPa', 'FontSize', 16)
yl = ylim; 
plot(med*[1 1], yl, 'r','LineWidth',2)
plot(pressure_thresh*[1 1], yl, 'k--','LineWidth',2)
xlim([0, 50])
title('Pressure Change 8/2018 -- 12/2020')

h2 = legend('Median',strcat( num2str(pressure_thresh), ' MPa'));
set(h2,'FontSize', 14)

% plot Survivor Function
[f,xi,bw] = ksdensity(dp,'support', 'positive','BoundaryCorrection',...
    'reflection', 'Function','Survivor');
display(['BandWidth:  ',num2str(bw)])
figure; area(xi,f,'facealpha',0.5,'edgecolor','none', 'HandleVisibility','off'); 
hold  on

set(gca,'FontSize', 14); hold on
xlabel('Pressure Change 8/2018 -- 12/2020, MPa', 'FontSize', 16)
ylabel('Probability', 'FontSize', 16)
xlim([0, 50])
yl = ylim; xl = xlim;
title('Survivor Function')

% to get interpolated probability for p > pressure_thresh
surv = interp1(xi,f,pressure_thresh);
plot([xl(1) pressure_thresh], surv*[1 1], 'r','LineWidth', 1.5)
plot(pressure_thresh*[1 1], [yl(1) surv], 'k--','LineWidth',2)

%% Account for the uncertainty in the threshold with uniform distribution
% over range [pressure_thresh - DeltaP_thresh, pressure_thresh + DeltaP_thresh]. 
% Choose Nthresh values over that range

DeltaP_thresh = 5; % in MPa
Nthresh = 11;
pressure_thresh_range = linspace(pressure_thresh - DeltaP_thresh, ...
    pressure_thresh + DeltaP_thresh, Nthresh);
DepthRange = 40*(pressure_thresh_range - pressure_thresh);

surv_range = zeros(Nthresh,1);
for k = 1:Nthresh
    surv_range(k) = interp1(xi,f,pressure_thresh_range(k));
end

 %plot results with range

figure
t = tiledlayout(1,1);
ax1 = axes(t);
plot(ax1,DepthRange,surv_range,'-b','LineWidth',2)

xlabel('Depth Relative to 2020 vent, m', 'FontSize', 16)
ylabel('Probability of Exceeding Threshold', 'FontSize', 16)
set(gca,'FontSize', 14); hold on
axis tight


ax2 = axes(t);
plot(ax2,pressure_thresh_range,surv_range,'LineWidth',2)
xlim = [pressure_thresh_range(1) pressure_thresh_range(end)];
ax2.XAxisLocation = 'top';
axis tight

ax2.YAxisLocation = 'left';
ax2.XTick = ax1.XTick/40 + 11.5;
xlabel('Pressure Threshold, MPa', 'FontSize', 16)
set(gca,'FontSize', 14); hold on
grid
yl = ylim;
plot(pressure_thresh*[1 1], yl, 'k:','LineWidth',2)

