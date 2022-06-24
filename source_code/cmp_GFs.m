function GFs = cmp_GFs(samples, gps_xy, mu, nu, V_HMM, V_SC, HMM_distr_flag)
% Compute Green's functions for sites with fixed HMM, SC locations, geometry,
% and orientations. 

% Input: 
% samples   = (NPDFSamples X Nparams) posterior samples of relevant
% gps_xy    = (3 X Nsites) measured GPS displacements
% mu        = shear modulus
% nu        = Poisson's ratio
% V_HMM     = volume of Halemaumau
% V_SC      = volume of South Caldera reservoir
% HMM_distr_flag = 1 for using Kyle's posterior distribution of HMM volume minus co-collapse volume reduction

% Output:
% GF        = cell of Green's functions 

Nsamples = length(samples(:,1));
dp_HMM = 1e6;
dp_SC = 1e6;

GFs = cell(Nsamples, 1);

for i = 1:Nsamples
    if HMM_distr_flag == 1
        V_HMM = samples(i, 11); % Take samples of HMM volume from Kyle's posterior distribution minus co-collapse volume reduction
    end
    
    dx_HMM = samples(i, 1); dy_HMM = samples(i, 2); dz_HMM = -samples(i, 3); alpha_HMM = samples(i, 4); 
    dx_SC = samples(i, 5); dy_SC = samples(i, 6); dz_SC = -samples(i, 7); alpha_SC = samples(i, 8); 
    dip_SC = samples(i, 9); 
    strike_SC = samples(i, 10); 

    a_HMM = (3*V_HMM/(4*pi)*alpha_HMM^2)^(1/3);
    b_HMM = a_HMM/alpha_HMM;

    a_SC = (3*V_SC/(4*pi)*alpha_SC^2)^(1/3);
    b_SC = a_SC/alpha_SC;

    % input for spheroid code
    m_HMM = [a_HMM; b_HMM; 90; 0; dx_HMM; dy_HMM; dz_HMM; -dp_HMM];
    m_SC = [a_SC; b_SC; dip_SC; strike_SC; dx_SC; dy_SC; dz_SC; -dp_SC];
    
    GF.HMM = spheroid(m_HMM,gps_xy,nu,mu);
    GF.SC = spheroid(m_SC,gps_xy,nu,mu);
    
    GFs{i} = GF;
end