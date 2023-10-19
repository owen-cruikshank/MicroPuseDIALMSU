function [absorption,cross_section] = absorption_O2_770_PCA_off(T,P,nu_Range,WV)
%File: absorption_O2_770_model_wavenumber.m
%Date: 02/28/2020
%Author: Owen Cruikshank
%Inputs:
%   -T:[K] scalar or (range x time) vector of atmopheric temperature as a function
%   of range
%   -P:[atm] scalar or (range x time) vector of atmopheric pressure as a function
%   of range
%   -nu_Range:[1/cm] wavenumber scan. Dimentions: (1 x 1 x wavenumber)
%   -WV:[1/m^3] water vapor number density. Dimentions: (range x time)
%
%Outputs:
%   -absorption: [1/m] the atmopheric absorption of O2 at the input
%   wavenumber arround the 780nm line, dimensions (range x time x wavenumber)
%   -sigma: [m^2] The voight absorption cross section of O2 at the input
%   wavenumber arround the 780nm line, dimensions (range x time x wavenumber)
%   -f: [m] Absorption lineshape function, dimensions (range x time x wavnumber)


%load variables needed

if length(nu_Range)==1
    load(fullfile('CalibrationData','PCA_1_11_21singleOnline.mat'),'M','muP','muT','muY','nu','sigmaP','sigmaT');
else
    %load(fullfile('CalibrationData','PCA_1_11_21single.mat'),'M','muP','muT','muY','nu','sigmaP','sigmaT');
    %load('PCA_8_2_22single.mat','M','muP','muT','muY','nu','sigmaP','sigmaT');
    %load('PCA_6_13_23single.mat','M','muP','muT','muY','nu','sigmaP','sigmaT');
    %load('PCA_6_19_23singleHITRAN2016.mat','M','muP','muT','muY','nu','sigmaP','sigmaT');

    %load('PCA_7_27_23singleHITRAN2016bigSpectrum.mat','M','muP','muT','muY','nu','sigmaP','sigmaT');
    load('PCA_7_31_23singleHITRAN2016bigSpectrum_off.mat','M','muP','muT','muY','nu','sigmaP','sigmaT');
end

%order
No = 20;
        normT = (T-muT)/sigmaT;
        normP = (P-muP)/sigmaP;
        
inc=1;
theta = ones(size(normT,1),size(normT,2),210);
for n = 1:No
    for m = 1:n
        theta(:,:,inc) = normT.^(n-m) .* normP.^(m-1);
        inc=inc+1;
    end
end

thetapermute = permute(theta(:,:,:),[3 2 1]);

cross_sectionI = ones(length(nu),length(T(1,:)),length(T(:,1)));
for j = 1:length(T(:,1))
    for i = 1:length(T(1,:))
        %cross_sectionI(:,i,j) = muY + M*thetapermute(:,i,j);  
        cross_sectionI(:,i,j) = M*thetapermute(:,i,j);  
    end
end
cross_sectionI = muY + cross_sectionI;  


cross_section = permute(cross_sectionI(:,:,:),[3 2 1]);

kB = 1.38065E-23;                       %[J/K] Boltzman's constant 
q_O2 = .2095;                           %[unitless] O2 atmospheric mixing ratio 
N_o2 = ((P*101325)./(kB*T)-WV) * q_O2;     %[molecule/m^3](t x r)O2 number density from atmopsheric number density and O2 mixing ratio

absorption = cross_section .* N_o2;     %[1/m](t x r x nu)absorption coefficeint of oxygen in the atmosphere
end