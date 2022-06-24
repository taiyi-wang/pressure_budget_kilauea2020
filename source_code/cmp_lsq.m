function [du, m] = cmp_lsq(GF, d)
% Least sqaures method to solve for the best fitting pressure change for
% the two reservoirs

% Input:
% GF = a structure with two fields: Halemaumau Green's functions and South
%      Caldera Green's functions
% d  = (Nsites X 3 , 1) data vector organized in the order of E, N, Z

% Output:
% du = (Nsites X 3 , 1) vector of best fit displacements in the order of E, N, Z
% m  = (2, 1) vector of HMM pressure change, SC pressure change

% 08/26/21

% reorganize the Green's functions into a matrix
G = [[GF.HMM(1,:)';GF.HMM(2,:)';GF.HMM(3,:)'],[GF.SC(1,:)';GF.SC(2,:)';GF.SC(3,:)']];

% Ignore entries of the Green's function and data vector corresponding to stations without data 
G_nonan = G(~isnan(d), :); 
d_nonan = d(~isnan(d));

% least squares
m = (G_nonan'*G_nonan)\G_nonan'*d_nonan;
m = m.*1e6; % recover pressure change in Pa

% predicted displacements from best fit model for ALL STATIONS
du = G*(m./1e6);
end