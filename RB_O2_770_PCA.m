function [sponts6] = RB_O2_770_PCA(T,P,nu_Range)
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
%Coherent s7
%load('RBPCA_1_15_21single.mat')
%spontaneous s6
%load(fullfile('CalibrationData','RBPCA_11_2_21single.mat'),'M','muP','muT','muY','nu','sigmaP','sigmaT')


%load('PCA_7_30_23singleRayleighBrillionbigSpectrum.mat','M','muP','muT','muY','nu','sigmaP','sigmaT')

%load('PCA_7_30_23singleRayleighBrillion.mat','M','muP','muT','muY','nu','sigmaP','sigmaT')


%load('PCA_10_06_23singleRayleighBrillion_offline.mat','muY','M','muT','muP','sigmaT','sigmaP','nu')
load(fullfile('CalibrationData','PCA_10_10_23singleRayleighBrillion_offline.mat'),'muY','M','muT','muP','sigmaT','sigmaP','nu')


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
sponts6I = ones(length(nu),length(T(1,:)),length(T(:,1)));
for j = 1:length(T(:,1))
    for i = 1:length(T(1,:))
        %sponts6I(:,i,j) = muY + M*thetapermute(:,i,j);  
        sponts6I(:,i,j) =  M*thetapermute(:,i,j);
    end
end
sponts6I = sponts6I +muY;

sponts6 = permute(sponts6I(:,:,:),[3 2 1]);
normsponts6 = trapz(sponts6,3)*(nu_Range(2)-nu_Range(1))*100;
sponts6 = sponts6./normsponts6;

end