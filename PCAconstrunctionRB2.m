%principal component analysis construction
function [Spectrum] = PCAconstrunctionRB2(Spectrum)
%addpath('C:\Users\Owen\OneDrive - Montana State University\Research\O2 DIAL\TempRetrievalRobert\PCASpectra\TestingCode')
N=200;
%N=1000;
%T vector
Tmin = 140;
Tmax = 320;
T = linspace(Tmin,Tmax,N);%[K]
%P vector
Pmin = .1;
Pmax = 1.2;
P = linspace(Pmin,Pmax,N);%[atm]

%freqency


%%
% lambda_online = 769.7958;
% lambda_offline = 770.1085;
% lambda_offline = Spectrum.lambda_offline(1);
% lambda_online = Spectrum.lambda_online(1);

% %lambda_online(1:length(ts)) = 769.7954;
% %lambda_offline(1:length(ts)) = 770.1081;
% nu_online = 10^7./lambda_online;                    %[cm-1] Online wavenumber
% nu_offline = 10^7./lambda_offline;                  %[cm-1] Offline wavenumber
% 
% nu01 = nu_online;                                   %[cm-1] Set center of scan to 
% nuMin = nu01-0.334;                                 %[cm-1] Scan lower bound
% nuMax = nu01+0.334;                                 %[cm-1] Scan upper bound
% nuBin = 0.00222;                                    %[cm-1] Scan increment
% %nuBin = 0.000222;                                    %[cm-1] Scan increment
% %nuBin = nuBin/2;
% %nuBin = 8.8800e-04;
% %nuBin = 0.0007; %[cm-1] Scan increment
% 
% % % % nuMin = nu_online-0.334*5;                                 %[cm-1] Scan lower bound
% % % % nuMax = nu_online+0.334*5;                                 %[cm-1] Scan upper bound
% % % % nuBin = 0.00222/3;  
%  nu_scan = (nuMin:nuBin:nuMax);                      %[cm-1](1 x nu) Scan vector
% 
% % nuMin = nu_online-0.334;                                 %[cm-1] Scan lower bound
% % nuMax = nu_online+0.334;                                 %[cm-1] Scan upper bound
% % nuBin = 0.00222;                                    %[cm-1] Scan increment
% % nu_scan = (nuMin:nuBin:nuMax);                      %[cm-1](1 x nu) Scan vector
% 
% nu01_off = nu_offline;
% nuMin_off = nu01_off-0.334;                                 %[cm-1] Scan lower bound
% nuMax_off = nu01_off+0.334;                                 %[cm-1] Scan upper bound
% 
% 
% % % % nuMin_off = nu01_off-0.334*5;                                 %[cm-1] Scan lower bound
% % % % nuMax_off = nu01_off+0.334*5;                                 %[cm-1] Scan upper bound
% nu_scan_off = (nuMin_off:nuBin:nuMax_off);
% 
% 
% 
% % %new nu
% % lambdaCenter = 769.7958;
% % lambdaWidth = 0.2;
% % lambdaDelta = .00002;
% % lambda = lambdaCenter-lambdaWidth/2:lambdaDelta:lambdaCenter+lambdaWidth/2;%[nm]
% % nu_scan = 10^7./lambda;
% % %nu_scan = permute(nu,[3 2 1]);
% 
% 
% %i_scan = length(nu_scan);                           %[none] length of scan vector
% 
% lambda_scan = 10^7./nu_scan;                        %[nm](1 x lambda) Scan vector
% lambda_scan_off = 10^7./nu_scan_off;
% %f_scan = nu_scan * c * 100;                         %[Hz]( x f) Scan vector
% 
% nu_scan_3D_short = permute(nu_scan, [3 1 2]);       %[cm-1] putting scan in third dimension
% nu_scan_3D_short_off = permute(nu_scan_off, [3 1 2]);       %[cm-1] putting scan in third dimension
% 
% lambda_scan_3D_short = 10^7./nu_scan_3D_short;
% lambda_scan_3D_short_off = 10^7./nu_scan_3D_short_off;
% i_scan_3D_short = length(nu_scan_3D_short);         %[none] length of scan vector
% 
% del_nu = nu_scan_3D_short-nu_online;                %[1/cm] difference from center
% del_lambda = lambda_scan_3D_short-lambda_online;
% nu = nu_scan_3D_short;
% lambda = lambda_scan;

% nu = permute(nu01,[3 2 1]);
% lambda = permute(lambda_online,[3 2 1]);
%%


% nu01 = 1.299045834338435e+04;                                   %[cm-1] Set center of scan to 
% nuMin = nu01-0.334;                                 %[cm-1] Scan lower bound
% nuMax = nu01+0.334;                                 %[cm-1] Scan upper bound
% nuBin = 0.00222;                                    %[cm-1] Scan increment
% %nuBin = 0.001;                                    %[cm-1] Scan increment
% nu_scan = (nuMin:nuBin:nuMax);                      %[cm-1](1 x nu) Scan vector
% %lambda = 10^7./nu_scan;
% %i_scan = length(nu_scan);                           %[none] length of scan vector
% 
% %nu = permute(nu_scan,[3 1 2]);
% nu = permute(nu,[3 1 2]);
% nuCenter = 10^7./lambdaCenter;


%construct P and T as j
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

%%
%create training matrix
%y = zeros(1,N*N,length(lambda));
% % % vw = zeros(N*N,1,1);
% % % q = 10;
% % % for i = 1:length(nu)
% % %     %for j = 1:N*N
% % %         [~,y(i,:),~] = absorption_O2_770_model(TP(:,1),TP(:,2),nu(i),0);
% % %         %[a,y,~] = absorption_O2_770_model_wavenumber(TP(:,1),TP(:,2),nu,0);
% % % 
% % %     %end
% % %     %y = permute(y,[3 2 1]);
% % % end

%%F
% %FreqSpec
% %FreqSpec = DLam2DNu(CenterLambda.*1e-9,(WavelengthsDesired-CenterLambda).*1e-9);
% %lambdaScan = permute(lambda_scan,[1 3 2]);
% FreqSpec = 299792458.*(lambda_scan-lambda_online).*1e-9./(lambda_online.*1e-9)./(lambda_online.*1e-9);
% for i = 1:length(TP(:,1))
%     parfor j = 1:length(TP(:,2))
% y(:,:,i,j) = RayleighBrillouinSpecWavelength(FreqSpec,lambda_online.*1e-9,TP(j,2).*101325,TP(i,1));
%     end
% end


% % FreqSpec = 299792458.*(lambda_scan-lambda_online).*1e-9./(lambda_online.*1e-9)./(lambda_online.*1e-9);
% % parfor i = 1:length(TP(:,1))
% %     y(i,:) = RayleighBrillouinSpecWavelength(FreqSpec,lambda_online.*1e-9,TP(i,2).*101325,TP(i,1));
% % end

lambda_offline = Spectrum.lambda_offline(1);
lambda_scan_off = permute(Spectrum.lambda_scan_3D_short_off,[3 2 1]);

FreqSpec = 299792458.*(lambda_scan_off-lambda_offline).*1e-9./(lambda_offline.*1e-9)./(lambda_offline.*1e-9);

TP1 = TP(:,1);
TP2 = TP(:,2);
parfor i = 1:length(TP(:,1))
    %y(i,:) = RayleighBrillouinSpecWavelength(FreqSpec,lambda_offline.*1e-9,TP2(i,1).*101325,TP1(i,1));
     y(i,:) = RayleighBrillouinSpecWavelength(FreqSpec,lambda_offline.*1e-9,TP(i,2).*101325,TP(i,1));
end

y=y';


%%
clear theta
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

%%
clear muY
muY = mean(y,2);

%%
clear U S V
[U,S,V] = svd(y);

%%
clear W
W = U'*(y-muY);

%%
clear C
C = W*pinv(theta);
%%
clear M
M = U*C;

Spectrum.RBoffline.muY =muY ;
Spectrum.RBoffline.M = M;
Spectrum.RBoffline.muT = muT;
Spectrum.RBoffline.muP = muP;
Spectrum.RBoffline.sigmaT = sigmaT;
Spectrum.RBoffline.sigmaP = sigmaP;



