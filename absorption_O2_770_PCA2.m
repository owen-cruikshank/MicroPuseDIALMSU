function [absorption,cross_section] = absorption_O2_770_PCA2(T,P,Spectrum,WV,Constants)
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


%==load variables needed
% if length(nu_Range)==1
%     load(fullfile('CalibrationData','PCA_1_11_21singleOnline.mat'),'M','muP','muT','muY','nu','sigmaP','sigmaT');
% else
%     %load(fullfile('CalibrationData','PCA_1_11_21single.mat'),'M','muP','muT','muY','nu','sigmaP','sigmaT');
%     %load('PCA_8_2_22single.mat','M','muP','muT','muY','nu','sigmaP','sigmaT');
%     %load('PCA_6_13_23single.mat','M','muP','muT','muY','nu','sigmaP','sigmaT');
%     load(fullfile('CalibrationData','PCA_6_19_23singleHITRAN2016.mat'),'M','muP','muT','muY','nu','sigmaP','sigmaT');
%     load(fullfile('CalibrationData','PCA_11_25_23singleHITRAN2012.mat'),'M','muP','muT','muY','nu','sigmaP','sigmaT');
% 
%     %load(fullfile('CalibrationData','PCA_10_30_23singleHITRAN2016doubleSpectrum.mat'),'M','muP','muT','muY','nu','sigmaP','sigmaT');
%     %load('PCA_7_27_23singleHITRAN2016bigSpectrum.mat','M','muP','muT','muY','nu','sigmaP','sigmaT');
% end

muY = Spectrum.Aonline.muY;
M = Spectrum.Aonline.M;
muT = Spectrum.Aonline.muT;
muP = Spectrum.Aonline.muP;
sigmaT = Spectrum.Aonline.sigmaT;
sigmaP = Spectrum.Aonline.sigmaP;

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

cross_sectionI = ones(length(Spectrum.nu_scan_3D_short),length(T(1,:)),length(T(:,1)));
for j = 1:length(T(:,1))
    for i = 1:length(T(1,:)) 
        cross_sectionI(:,i,j) = M*thetapermute(:,i,j);  
    end
end
cross_sectionI = muY + cross_sectionI;  


cross_section = permute(cross_sectionI(:,:,:),[3 2 1]);

N_o2 = ((P*Constants.ATMtoPA)./(Constants.kB*T)-WV) * Constants.q_O2;
absorption = cross_section .* N_o2;     %[1/m](t x r x nu)absorption coefficeint of oxygen in the atmosphere
end