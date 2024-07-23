%principal component analysis construction
function [Spectrum] = PCAconstrunctionRB2(Spectrum,Constant)

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

%Calulate offline
lambda_offline = Spectrum.lambda_offline(1);
lambda_scan_off = permute(Spectrum.lambda_scan_3D_short_off,[3 2 1]);
%Convert lambda to freqency vector
FreqSpecOffline = 299792458.*(lambda_scan_off-lambda_offline).*1e-9./(lambda_offline.*1e-9)./(lambda_offline.*1e-9);

lambda_online = Spectrum.lambda_online(1);
lambda_scan = permute(Spectrum.lambda_scan_3D_short,[3 2 1]);
%Convert lambda to freqency vector
FreqSpecOnline = 299792458.*(lambda_scan-lambda_online).*1e-9./(lambda_online.*1e-9)./(lambda_online.*1e-9);

if FreqSpecOnline(Spectrum.online_index)==0
FreqSpecOnline(Spectrum.online_index) = (FreqSpecOnline(Spectrum.online_index)-FreqSpecOnline(Spectrum.online_index-1))/10;
FreqSpecOffline(Spectrum.offline_index) = (FreqSpecOffline(Spectrum.offline_index)-FreqSpecOffline(Spectrum.offline_index-1))/10;
end

yRBOff = zeros(length(TP(:,1)),length(FreqSpecOffline));
yRBOn = zeros(length(TP(:,1)),length(FreqSpecOnline));
%Calculate training Spectrum
parfor i = 1:length(TP(:,1))
     yRBOff(i,:) = RayleighBrillouinSpecWavelength(FreqSpecOnline,lambda_offline.*1e-9,TP(i,2).*101325,TP(i,1));
     yRBOn(i,:) = RayleighBrillouinSpecWavelength(FreqSpecOffline,lambda_online.*1e-9,TP(i,2).*101325,TP(i,1));
     
end

yAOn=zeros(length(TP(:,1)),1,length(FreqSpecOnline));
for i = 1:length(Spectrum.nu_scan_3D_short)
    [~,yAOn(:,:,i)] = absorption_O2_770_model(TP(:,1),TP(:,2),Spectrum.nu_scan_3D_short(i),0,Constant);
end
yAOn=permute(yAOn,[1 3 2]);
yRBOff=yRBOff';
yRBOn=yRBOn';
yAOn=yAOn';


%Create theta based on T and P
theta = zeros(210,N*N);
No = 20;
for j = 1:N*N
    normT = (TP(j,1)-muT)/sigmaT;
    normP = (TP(j,2)-muP)/sigmaP;
    inc=1;
    for n = 1:No
        for m = 1:n
            theta(inc,j) = normT^(n-m) * normP^(m-1);
            inc=inc+1;
        end
    end
end

% == Create PCA variables Rayleigh brillion offline
muY = mean(yRBOff,2);
[U,S,V] = svd(yRBOff); %Varibles S and V MUST be specified. I am not sure why
clear S V
W = U'*(yRBOff-muY);
C = W*pinv(theta);
M = U*C;
%Variables needed for PCA reconstruction
Spectrum.RBoffline.muY =muY ;
Spectrum.RBoffline.M = M;
Spectrum.RBoffline.muT = muT;
Spectrum.RBoffline.muP = muP;
Spectrum.RBoffline.sigmaT = sigmaT;
Spectrum.RBoffline.sigmaP = sigmaP;
%%
% == Create PCA variables Rayleigh brillion online
muY = mean(yRBOn,2);
[U,S,V] = svd(yRBOn); %Varibles S and V MUST be specified. I am not sure why
clear S V
W = U'*(yRBOn-muY);
C = W*pinv(theta);
M = U*C;
Spectrum.RBonline.muY =muY ;
Spectrum.RBonline.M = M;
Spectrum.RBonline.muT = muT;
Spectrum.RBonline.muP = muP;
Spectrum.RBonline.sigmaT = sigmaT;
Spectrum.RBonline.sigmaP = sigmaP;
%%
%== Create PCA variables for online absorption
muY = mean(yAOn,2);
[U,S,V] = svd(yAOn); %Varibles S and V MUST be specified. I am not sure why
clear S V
W = U'*(yAOn-muY);
C = W*pinv(theta);
M = U*C;
Spectrum.Aonline.muY =muY ;
Spectrum.Aonline.M = M;
Spectrum.Aonline.muT = muT;
Spectrum.Aonline.muP = muP;
Spectrum.Aonline.sigmaT = sigmaT;
Spectrum.Aonline.sigmaP = sigmaP;



