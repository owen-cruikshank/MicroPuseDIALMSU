function [sponts6] = RB_O2_770_PCA(T,P,nu_Range,Spectrum)
%File: absorption_O2_770_model_wavenumber.m
%Date: 02/28/2020
%Author: Owen Cruikshank
%Inputs:
%   -T:[K] scalar or (range x time) vector of atmopheric temperature as a function
%   of range
%   -P:[atm] scalar or (range x time) vector of atmopheric pressure as a function
%   of range
%   -nu_Range:[1/cm] wavenumber scan. Dimentions: (1 x 1 x wavenumber)
%
%Outputs:


muY = Spectrum.RBoffline.muY;
M = Spectrum.RBoffline.M;
muT = Spectrum.RBoffline.muT;
muP = Spectrum.RBoffline.muP;
sigmaT = Spectrum.RBoffline.sigmaT;
sigmaP = Spectrum.RBoffline.sigmaP;

%Recaculate theta with input T and P
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
sponts6I = ones(length(Spectrum.nu_scan_3D_short_off),length(T(1,:)),length(T(:,1)));
for j = 1:length(T(:,1))
    for i = 1:length(T(1,:))
        sponts6I(:,i,j) = muY + M*thetapermute(:,i,j);  
    end
end

%Move spectrum to third dimension
sponts6 = permute(sponts6I(:,:,:),[3 2 1]);
%Normalize spectrum to integral
normsponts6 = trapz(sponts6,3)*(nu_Range(2)-nu_Range(1))*100;
sponts6 = sponts6./normsponts6;

end