%principal component analysis construction
function [Spectrum] = PCAconstrunctionRB2(Spectrum)

% Input varibles
%Spectrum.lambda_offline(1); -> Offline wavelength [nm]
%Spectrum.lambda_scan_3D_short_off -> vector for offline evaluation spectrum [nm] along third dimension.

%Creating Temperature and Pressure vectors for evaulation
N=200;
%T vector
Tmin = 140;
Tmax = 320;
T = linspace(Tmin,Tmax,N);%[K]
%P vector
Pmin = .1;
Pmax = 1.2;
P = linspace(Pmin,Pmax,N);%[atm]


%Construct P and T as j vector
TP = zeros(N*N,2);
increment = 1;
for i = 1:N
    for j = 1:N
        TP(increment,1) = T(i);
        TP(increment,2) = P(j);
        increment = increment+1;
    end
end
muT = mean(T);
sigmaT = std(T);
muP = mean(P);
sigmaP = std(P);


lambda_offline = Spectrum.lambda_offline(1);
lambda_scan_off = permute(Spectrum.lambda_scan_3D_short_off,[3 2 1]);
%Convert lambda to freqency vector
FreqSpec = 299792458.*(lambda_scan_off-lambda_offline).*1e-9./(lambda_offline.*1e-9)./(lambda_offline.*1e-9);

%Calculate training Spectrum
parfor i = 1:length(TP(:,1))
     y(i,:) = RayleighBrillouinSpecWavelength(FreqSpec,lambda_offline.*1e-9,TP(i,2).*101325,TP(i,1));
end

y=y';


%Create theta based on T and P
No = 20;
for j = 1:N*N
    normT = (TP(j,1)-muT)/sigmaT;
    normP = (TP(j,2)-muP)/sigmaP;
    inc=1;
    for n = 1:No
        for m = 1:n
            theta(j,inc) = normT^(n-m) * normP^(m-1);
            inc=inc+1;
        end
    end
end
theta = theta';

muY = mean(y,2);

[U,S,V] = svd(y); %Varibles S and V MUST be specified. I am not sure why
clear S V

W = U'*(y-muY);
C = W*pinv(theta);
M = U*C;

%Variables needed for PCA reconstruction

Spectrum.RBoffline.muY =muY ;
Spectrum.RBoffline.M = M;
Spectrum.RBoffline.muT = muT;
Spectrum.RBoffline.muP = muP;
Spectrum.RBoffline.sigmaT = sigmaT;
Spectrum.RBoffline.sigmaP = sigmaP;



