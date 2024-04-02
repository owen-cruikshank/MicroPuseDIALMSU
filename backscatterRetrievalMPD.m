function [HSRL] = backscatterRetrievalMPD(Counts, Model, Spectrum, Options)
%File: backscatterRetrievalMPD.m
%Author: Owen Cruikshank
%Inputs:
%   -Counts: structure for measured counts
%   -Model: strucure for model T, P, and WV
%   -Spectrum: structure for spectrum used in calculation
%   -Options:structure for options
%
%Outputs:
%   -HSRL: structure for HSRL
%
%Description:
%HSRL retrieval for the MPD


[sponts6_off] = RB_O2_770_PCA(Model.T,Model.P,Spectrum.nu_scan_3D_short_off,Spectrum);
sponts6_off = sponts6_off./trapz(permute(Spectrum.nu_scan_3D_short_off,[3 2 1]),sponts6_off,3);

if strcmp( Options.MPDname,'03')
    load(fullfile('CalibrationData','MPD03ScanData.mat'),'CalInfo')
    CalInfo.ScanData.O2OfflineMol.Wavelength = CalInfo.ScanData.O2OfflineMol.Wavelength+Spectrum.WavemeterOffset;
    CalInfo.ScanData.O2OfflineComb.Wavelength = CalInfo.ScanData.O2OfflineComb.Wavelength+Spectrum.WavemeterOffset;
    CalInfo.ScanData.O2OnlineMol.Wavelength = CalInfo.ScanData.O2OnlineMol.Wavelength+Spectrum.WavemeterOffset;
    CalInfo.ScanData.O2OnlineComb.Wavelength = CalInfo.ScanData.O2OnlineComb.Wavelength+Spectrum.WavemeterOffset;

    CalInfo.ScanData.WVOffline.Wavelength = CalInfo.ScanData.WVOffline.Wavelength+Spectrum.WavemeterOffset;
    CalInfo.ScanData.WVOnline.Wavelength = CalInfo.ScanData.WVOnline.Wavelength+Spectrum.WavemeterOffset;
    
    offlineMolecularTransmission = interp1(10^7./CalInfo.ScanData.O2OfflineMol.Wavelength,  CalInfo.ScanData.O2OfflineMol.Transmission,  Spectrum.nu_scan_3D_short_off);
    offlineCombinedTransmission  = interp1(10^7./CalInfo.ScanData.O2OfflineComb.Wavelength, CalInfo.ScanData.O2OfflineComb.Transmission, Spectrum.nu_scan_3D_short_off);
    onlineMolecularTransmission  = interp1(10^7./CalInfo.ScanData.O2OnlineMol.Wavelength,   CalInfo.ScanData.O2OnlineMol.Transmission,   Spectrum.nu_scan_3D_short);
    onlineCombinedTransmission   = interp1(10^7./CalInfo.ScanData.O2OnlineComb.Wavelength,  CalInfo.ScanData.O2OnlineComb.Transmission,  Spectrum.nu_scan_3D_short);

    WVOnlineTransmission = interp1(10^7./CalInfo.ScanData.WVOnline.Wavelength,  CalInfo.ScanData.WVOnline.Transmission,  Spectrum.nu_scanwv_3D_short);
elseif strcmp(Options.MPDname,'00')
    load('CalibrationData\TransmissionData20220809.mat','Data_Wavelength')

    onlineCombinedTransmission   = interp1((Data_Wavelength.lambda_on.*10^9) +Spectrum.WavemeterOffset,Data_Wavelength.Nc_on, Spectrum.lambda_scan_3D_short);
    onlineMolecularTransmission  = interp1((Data_Wavelength.lambda_on.*10^9) +Spectrum.WavemeterOffset,Data_Wavelength.Nm_on, Spectrum.lambda_scan_3D_short);
    offlineCombinedTransmission  = interp1((Data_Wavelength.lambda_off.*10^9)+Spectrum.WavemeterOffset,Data_Wavelength.Nc_off,Spectrum.lambda_scan_3D_short_off);
    offlineMolecularTransmission = interp1((Data_Wavelength.lambda_off.*10^9)+Spectrum.WavemeterOffset,Data_Wavelength.Nm_off,Spectrum.lambda_scan_3D_short_off);
    
    WVOnlineTransmission = ones(size(offlineMolecularTransmission));

elseif strcmp(Options.MPDname,'05')
    load('CalibrationData\CalibrationTablesBoulderSponS6062021.mat','BoulderHSRLcoefficentsSponS6_062021')

    BoulderHSRLcoefficents062021 = BoulderHSRLcoefficentsSponS6_062021;
    P=BoulderHSRLcoefficents062021.P; 
    T=permute(BoulderHSRLcoefficents062021.T,[2 1]);
    Eta_m=1;
    Eta_c=1;
    Cmm=permute(BoulderHSRLcoefficents062021.Cmm,[3 2 1]);
    Cmc=permute(BoulderHSRLcoefficents062021.Cmc,[3 2 1]); 
    Cam=BoulderHSRLcoefficents062021.Cam; 
    Cac=BoulderHSRLcoefficents062021.Cac;
    
    
    Cmm2=zeros(size(Counts.o2on));
    Cmc2=zeros(size(Counts.o2on));

    wdP = Model.P;
    wdT = Model.T;

%     for i=1:length(LidarData.Time)
%        parfor j=1:length(LidarData.Range)
%         Cmm2(j,i)=interp2(P,T,Cmm,wdP(j,i),wdT(j,i));
%         Cmc2(j,i)=interp2(P,T,Cmc,wdP(j,i),wdT(j,i));
%        end
%     end

    for j=1:size(Counts.o2on,1)
    Cmm2(j,:)=interp2(P,T,Cmm,wdP(j,:),wdT(j,:));
    Cmc2(j,:)=interp2(P,T,Cmc,wdP(j,:),wdT(j,:));
    end

    Cmm = Cmm2;
    Cmc = Cmc2;
    
    %Calibration Constants
    LidarData.Cac=Cac; %Aerosol in Combined
    LidarData.Cam=Cam; %Aerosol in Molecular

    %Beamsplitter and Detector Efficiencies
    etac=Eta_c;
    etam=Eta_m;
        eta_c=Eta_c;
    eta_m=Eta_m;

    onlineCombinedTransmission   = ones(size(Spectrum.lambda_scan_3D_short));
    onlineMolecularTransmission  = ones(size(Spectrum.lambda_scan_3D_short));
    offlineCombinedTransmission  = ones(size(Spectrum.lambda_scan_3D_short));
    offlineMolecularTransmission = ones(size(Spectrum.lambda_scan_3D_short));
    
    WVOnlineTransmission = ones(size(Spectrum.nu_scanwv_3D_short));
end

HSRL.offlineMolecularTransmission = offlineMolecularTransmission;
HSRL.offlineCombinedTransmission = offlineCombinedTransmission;
HSRL.onlineMolecularTransmission = onlineMolecularTransmission;
HSRL.onlineCombinedTransmission = onlineCombinedTransmission;
HSRL.WVOnlineTransmission =WVOnlineTransmission;

if strcmp(Options.MPDname,'00') || strcmp(Options.MPDname,'03')
int = [1:50 size(offlineMolecularTransmission,3)-50:size(offlineMolecularTransmission,3)];

eta_m=(mean(offlineCombinedTransmission(:,:,int)./offlineMolecularTransmission(:,:,int))+1).^-1;

% etam = (mean(CalInfo.ScanData.O2OnlineComb.Transmission./CalInfo.ScanData.O2OnlineMol.Transmission)+1)^-1;
 %eta_m = (mean(onlineCombinedTransmission./onlineMolecularTransmission)+1)^-1;
 eta_c = 1-eta_m;

 etam = eta_m;
 etac = eta_c;


% %normalize

offlineMolecularTransmission = offlineMolecularTransmission./(max(offlineCombinedTransmission) .*etam./etac);
offlineCombinedTransmission  = offlineCombinedTransmission ./max(offlineCombinedTransmission);

aerosolSpectrum = zeros(size(Spectrum.nu_scan_3D_short_off));
aerosolSpectrum(:,:,Spectrum.offline_index(1))=1;
aerosolSpectrum = aerosolSpectrum./trapz(permute(Spectrum.nu_scan_3D_short_off,[3 2 1]),aerosolSpectrum,3);

Cmm = trapz(permute(Spectrum.nu_scan_3D_short_off,[3 2 1]),sponts6_off.*offlineMolecularTransmission,3);
Cmc = trapz(permute(Spectrum.nu_scan_3D_short_off,[3 2 1]),sponts6_off.*offlineCombinedTransmission,3);
Cam = trapz(permute(Spectrum.nu_scan_3D_short_off,[3 2 1]),aerosolSpectrum.*offlineMolecularTransmission,3);
Cac = trapz(permute(Spectrum.nu_scan_3D_short_off,[3 2 1]),aerosolSpectrum.*offlineCombinedTransmission,3);
end

Cmm = Cmm./Cac;
Cam = Cam./Cac;
Cmc = Cmc./Cac;
Cac = Cac./Cac;

%Bm = 5.45e-32 .* Model.P.*101325./Model.T./1.380649e-23.*550^4./(Spectrum.lambda_offline(1).^4);

Bm = 9.94266e-7*(Model.P./0.009869233)   ./(Model.T)   ;

%Gm=1;
% Sm = (Gm.*Counts.o2off_mol-Cam.*Counts.o2off)./(Cmm-Cam);
% Ba = (Counts.o2off./Sm-1).*Bm;
% 


%Sm = (Gm.*Counts.Nm_off-Cam.*Counts.Nc_off)./(Cmm-Cam);
%Ba = (Counts.Nc_off./Sm-1).*Bm;

BSR = 1 - (Counts.Nm_on.*Cmm.*Counts.Nc_off-Counts.Nc_on.*Cmc.*Counts.Nm_off)./(Counts.Nm_on.*Cam.*Counts.Nc_off-Counts.Nc_on.*Counts.Nm_off);


Ba = Bm .* (BSR-1);

%Output variables
HSRL.BSR=BSR;
HSRL.Ba = Ba;
HSRL.Bm = Bm;
HSRL.eta_m = eta_m;
HSRL.eta_c = eta_c;
HSRL.Cam = Cam;
HSRL.Cmm = Cmm;
HSRL.Cac = Cac;
HSRL.Cmc = Cmc;


