function [HSRL] = backscatterRetrievalMPD03(Counts, Model, Spectrum, Options)
% % Model.T= 240;
% % Model.P = 71*0.00986923;

%[sponts6] = RB_O2_770_PCA(Model.T,Model.P,Spectrum.nu_scan_3D_short);
[sponts6_off] = RB_O2_770_PCA(Model.T,Model.P,Spectrum.nu_scan_3D_short_off,Spectrum);

%[sponts6_off] = RB_O2_770_PCA(Model.Tlong,Model.Plong,Spectrum.nu_scan_3D_short_off);
sponts6_off = sponts6_off./trapz(permute(Spectrum.nu_scan_3D_short_off,[3 2 1]),sponts6_off,3);

if strcmp( Options.MPDname,'03')
    %load('C:\Users\Owen\OneDrive - Montana State University\Research\O2 DIAL\Data\NCAR Boulder Data\MPD03ScanData.mat')
    load(fullfile('CalibrationData','MPD03ScanData.mat'),'CalInfo')
     CalInfo.ScanData.O2OfflineMol.Wavelength = CalInfo.ScanData.O2OfflineMol.Wavelength+Spectrum.WavemeterOffset;
     CalInfo.ScanData.O2OfflineComb.Wavelength = CalInfo.ScanData.O2OfflineComb.Wavelength+Spectrum.WavemeterOffset;
      CalInfo.ScanData.O2OnlineMol.Wavelength = CalInfo.ScanData.O2OnlineMol.Wavelength+Spectrum.WavemeterOffset;
     CalInfo.ScanData.O2OnlineComb.Wavelength = CalInfo.ScanData.O2OnlineComb.Wavelength+Spectrum.WavemeterOffset;
    
     %  CalInfo.ScanData.O2OfflineMol.Wavelength = CalInfo.ScanData.O2OfflineMol.Wavelength-Spectrum.WavemeterOffset;
     % CalInfo.ScanData.O2OfflineComb.Wavelength = CalInfo.ScanData.O2OfflineComb.Wavelength-Spectrum.WavemeterOffset;
     %  CalInfo.ScanData.O2OnlineMol.Wavelength = CalInfo.ScanData.O2OnlineMol.Wavelength-Spectrum.WavemeterOffset;
     % CalInfo.ScanData.O2OnlineComb.Wavelength = CalInfo.ScanData.O2OnlineComb.Wavelength-Spectrum.WavemeterOffset;
    
    offlineMolecularTransmission = interp1(10^7./CalInfo.ScanData.O2OfflineMol.Wavelength, CalInfo.ScanData.O2OfflineMol.Transmission, Spectrum.nu_scan_3D_short_off);
    offlineCombinedTransmission = interp1(10^7./CalInfo.ScanData.O2OfflineComb.Wavelength, CalInfo.ScanData.O2OfflineComb.Transmission, Spectrum.nu_scan_3D_short_off);
    onlineMolecularTransmission = interp1(10^7./CalInfo.ScanData.O2OnlineMol.Wavelength, CalInfo.ScanData.O2OnlineMol.Transmission, Spectrum.nu_scan_3D_short);
    onlineCombinedTransmission = interp1(10^7./CalInfo.ScanData.O2OnlineComb.Wavelength, CalInfo.ScanData.O2OnlineComb.Transmission, Spectrum.nu_scan_3D_short);

elseif strcmp(Options.MPDname,'00')
    load('CalibrationData\TransmissionData20220809.mat','Data_Wavelength')

onlineCombinedTransmission = interp1((Data_Wavelength.lambda_on.*10^9)+Spectrum.WavemeterOffset,Data_Wavelength.Nc_on,Spectrum.lambda_scan_3D_short);
onlineMolecularTransmission = interp1((Data_Wavelength.lambda_on.*10^9)+Spectrum.WavemeterOffset,Data_Wavelength.Nm_on,Spectrum.lambda_scan_3D_short);
offlineCombinedTransmission = interp1((Data_Wavelength.lambda_off.*10^9)+Spectrum.WavemeterOffset,Data_Wavelength.Nc_off,Spectrum.lambda_scan_3D_short_off);
offlineMolecularTransmission = interp1((Data_Wavelength.lambda_off.*10^9)+Spectrum.WavemeterOffset,Data_Wavelength.Nm_off,Spectrum.lambda_scan_3D_short_off);

end

HSRL.offlineMolecularTransmission = offlineMolecularTransmission;
HSRL.offlineCombinedTransmission = offlineCombinedTransmission;
HSRL.onlineMolecularTransmission = onlineMolecularTransmission;
HSRL.onlineCombinedTransmission = onlineCombinedTransmission;

% etam = (mean(CalInfo.ScanData.O2OnlineComb.Transmission./CalInfo.ScanData.O2OnlineMol.Transmission)+1)^-1;
 etam = (mean(onlineCombinedTransmission./onlineMolecularTransmission)+1)^-1;
 etac = 1-etam;
% 
% 
% %normalize
% etam = 0.5649;
% etac = 0.4351;
offlineMolecularTransmission = offlineMolecularTransmission./(max(offlineCombinedTransmission) .*etam./etac);
offlineCombinedTransmission = offlineCombinedTransmission./max(offlineCombinedTransmission);
% 





aerosolSpectrum = zeros(size(Spectrum.nu_scan_3D_short_off));
aerosolSpectrum(:,:,Spectrum.offline_index(1))=1;
aerosolSpectrum = aerosolSpectrum./trapz(permute(Spectrum.nu_scan_3D_short_off,[3 2 1]),aerosolSpectrum,3);

Cmm = trapz(permute(Spectrum.nu_scan_3D_short_off,[3 2 1]),sponts6_off.*offlineMolecularTransmission,3);
Cmc = trapz(permute(Spectrum.nu_scan_3D_short_off,[3 2 1]),sponts6_off.*offlineCombinedTransmission,3);
Cam = trapz(permute(Spectrum.nu_scan_3D_short_off,[3 2 1]),aerosolSpectrum.*offlineMolecularTransmission,3);
Cac = trapz(permute(Spectrum.nu_scan_3D_short_off,[3 2 1]),aerosolSpectrum.*offlineCombinedTransmission,3);


Cmm = Cmm./Cac;
Cam = Cam./Cac;
Cmc = Cmc./Cac;
%Cac = Cac./Cac;

%Bm = 5.45e-32 .* Model.P.*101325./Model.T./1.380649e-23.*550^4./(Spectrum.lambda_offline(1).^4);

Bm = 9.94266e-7*(Model.P./0.009869233)   ./(Model.T)   ;
%Bm = 9.94266e-7*(Model.Plong./0.009869233)   ./(Model.Tlong)   ;

%Gm=1;
% Sm = (Gm.*Counts.o2off_mol-Cam.*Counts.o2off)./(Cmm-Cam);
% Ba = (Counts.o2off./Sm-1).*Bm;
% 
% BSR2 = 1- (Counts.o2on_mol.*Cmm.*Counts.o2off-Counts.o2on.*Cmc.*Counts.o2off_mol)./(Counts.o2on_mol.*Cam.*Counts.o2off-Counts.o2on.*Counts.o2off_mol);


%Sm = (Gm.*Counts.Nm_off-Cam.*Counts.Nc_off)./(Cmm-Cam);
%Ba = (Counts.Nc_off./Sm-1).*Bm;

%BSR2 = 1- (Counts.Nm_on.*Cmm.*Counts.Nc_off-Counts.o2on.*Cmc.*Counts.Nm_off)./(Counts.Nm_on.*Cam.*Counts.Nc_off-Counts.Nc_on.*Counts.Nm_off);
BSR2 = 1 - (Counts.Nm_on.*Cmm.*Counts.Nc_off-Counts.Nc_on.*Cmc.*Counts.Nm_off)./(Counts.Nm_on.*Cam.*Counts.Nc_off-Counts.Nc_on.*Counts.Nm_off);

%BSR = (Ba+Bm)./Bm;

BSR = BSR2;

Ba = Bm .* (BSR-1);

HSRL.BSR=BSR;
HSRL.Ba = Ba;
HSRL.Bm = Bm;

5;

