function [cross_section] = cross_section_wv_828_PCA(T,P,~)
%File: absorption_O2_770_model_wavenumber.m
%Date: 03/22/2022
%Author: Owen Cruikshank
%Inputs:
%   -T:[K] scalar or (range x time) vector of atmopheric temperature as a function
%   of range
%   -P:[atm] scalar or (range x time) vector of atmopheric pressure as a function
%   of range
%   -nu_Range:[1/cm] wavenumber scan. Dimentions: (1 x 1 x wavenumber)
%
%Outputs:
%   -sigma: [m^2] The voight absorption cross section of O2 at the input
%   wavenumber arround the 780nm line, dimensions (range x time x wavenumber)

%load(fullfile('CalibrationData','PCA_3_15_22WV828single.mat'),'M','muP','muT','muY','nu','sigmaP','sigmaT')

load(fullfile('CalibrationData','PCA_8_10_22WV828single.mat'),'M','muP','muT','muY','nu','sigmaP','sigmaT')

%order
No = 20;
        normT = (T-muT)/sigmaT;
        normP = (P-muP)/sigmaP;
        
inc=1;
theta = ones(size(normT,1),size(normT,2),210);
for n = 1:No
    for m = 1:n
        theta(:,:,inc) = normT(:,:).^(n-m) .* normP(:,:).^(m-1);
        inc=inc+1;
    end
end

thetapermute = permute(theta(:,:,:),[3 2 1]);

cross_sectionI = ones(length(nu),length(T(1,:)),length(T(:,1)));
for j = 1:length(T(:,1))
    for i = 1:length(T(1,:))
        cross_sectionI(:,i,j) = muY + M*thetapermute(:,i,j);  
    end
end

cross_section = permute(cross_sectionI(:,:,:),[3 2 1]);

end