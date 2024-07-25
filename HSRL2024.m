%2024 HSRL
%7_24_24
%Owen Cruikshank
clear all
%====================
%==== Constants =====
%====================
Constants.g0 = 9.80665;       %[m/s^2] Gravitational acceleration 
Constants.M_air = 0.0289644;  %[kg/mol] Molar mass of Earth's air 
Constants.R = 8.3144598;      %[J/(mol*K)] Universal gas constant 

Constants.c = 2.99792458E8;           %[m/s] Speed of light 
Constants.kB = 1.38065E-23;           %[J/K][m^2 kg s-2 K-1] Boltzman's constant 
Constants.h = 6.62607004E-34;              %[Js] Planck's constant 
Constants.mo2 = 5.314E-26;            %[kg] Mass O2 molecule 
Constants.mWV = 2.991577548987048e-26;           %[kg] Mass H2O molecule
Constants.m_air = 4.792E-26;          %[kg] Mass of air
Constants.q_O2 = .2095;               %[unitless] O2 atmospheric mixing ratio
Constants.No = 2.47937E25;            %[1/m^3] Loschmidt's number  (referenced to 296 K and 1 atm)
Constants.CtoK = 273.15;             %Converstion from Celcius to Kelvin
Constants.mbtoAtm = 1/1013.25;         %converstion from milibars to atmopheres 
Constants.ATMtoPA = 101325;          %Converstion from atm to pa [kg m^-1 s^-2/atm]



%load data counts
spanDays = datetime(2024,7,22,'TimeZone','UTC');%yyyy,mm,dd
SpanDays1 = datetime(2024,7,23,'TimeZone','UTC');%yyyy,mm,dd
SpanDays2 = datetime(2024,7,23,'TimeZone','UTC');%yyyy,mm,dd
spanDays = SpanDays1:SpanDays2;
Options.BinTotal = 490;
Options.intTime = 30;
Options.path = fullfile('E:','Data','NetCDFOutput'); %Path for instument data
Options.MPDname = '000';
Data = loadMSUNETcdf(spanDays,Options);

Time = Data.MCS.Channel0.NewTimeGrid;

Range.nsPerBin = 250; %[ns] bin length in nanosections
Range.NBins = floor(Options.BinTotal); %number of range bins in vector
Range.rangeBin = (Constants.c * Range.nsPerBin(1)*10^-9)/2; %(m)range bin length
%Createraw range vector with offset for pulse length
Range.rm_raw_o2 = -Range.rangeBin*2:Range.rangeBin:Range.NBins(1)*Range.rangeBin-Range.rangeBin*3;
Range.rm_raw_o2 = Range.rm_raw_o2(:);                           %[m] Convert range vector to column vector

%load Transmission

data = ncread(fullfile(Options.path,'20240723','ReceiverScanMCS_00_20240723_190000.nc'),'Data');
datatime = ncread(fullfile(Options.path,'20240723','ReceiverScanMCS_00_20240723_190000.nc'),'time');
wavelengthtime = ncread(fullfile(Options.path,'20240723','ReceiverScanWavemeter_00_20240723_190000.nc'),'time');
wavelength = ncread(fullfile(Options.path,'20240723','ReceiverScanWavemeter_00_20240723_190000.nc'),'Wavelength');

wavelengthInterp = interp1(wavelengthtime,wavelength,datatime);
dataMol = mean(data(1:4:end,:),2);
dataComb = mean(data(2:4:end,:),2);

[wavelengthInterp,I] = sort(wavelengthInterp(1:4:end));
dataMol = movmean(dataMol(I),10);
dataComb = movmean(dataComb(I),10);

[wavelengthInterp,IA] = unique(wavelengthInterp);
dataMol = dataMol(IA);
dataComb = dataComb(IA);
wavelengthInterp = wavelengthInterp(~isnan(wavelengthInterp));
dataMol = dataMol(~isnan(wavelengthInterp));
dataComb = dataComb(~isnan(wavelengthInterp));


%load surface temperature and pressure

Press = 1;
Temp = 296;

Tr = Temp+ Range.rm_raw_o2.*(-6.5/1000);
gamma = Constants.g0 * Constants.M_air / Constants.R;     %[K/m]gravity molar mass of air and gas constant
Pr = Press.*exp(-cumtrapz(Range.rm_raw_o2,gamma./Tr,1));

%create spectrum
lambdaCentral = 780.2464;
nuCentralon = 10^7./lambdaCentral;
nuMin = nuCentralon-0.334;                                 %[cm-1] Scan lower bound
nuMax = nuCentralon+0.334;                                 %[cm-1] Scan upper bound
Spectrum.nuBin = 0.00222;                                    %[cm-1] Scan increment
Spectrum.nuBin = 0.00222/2;   
nu_scan = (nuMin:Spectrum.nuBin:nuMax);                      %[cm-1](1 x nu) Scan vector
Spectrum.nu_scan_3D_short = permute(nu_scan, [3 1 2]);       %[cm-1] putting scan in third dimension
lambda_scan = 10^7./nu_scan;
[~,Spectrum.offline_index] = min(abs(Spectrum.nu_scan_3D_short-nuCentralon));
%create RB spectrum
FreqSpec = 299792458.*(lambda_scan-lambdaCentral).*1e-9./(lambdaCentral.*1e-9)./(lambdaCentral.*1e-9);
for i = 1:length(Range.rm_raw_o2)

    [sponts6_off(i,:,:)] = RayleighBrillouinSpecWavelength(FreqSpec,lambdaCentral.*1e-9,Pr(i,:)*101325,Tr(i,:));
end
%sponts6_off = permute(sponts6_off, [3 1 2]);

sponts6_off = sponts6_off./trapz(permute(Spectrum.nu_scan_3D_short,[3 2 1]),sponts6_off,3);
%create HSRL calibration coefficients

offlineMolecularTransmission = interp1(wavelengthInterp,dataMol,lambda_scan);
offlineCombinedTransmission = interp1(wavelengthInterp,dataComb,lambda_scan);
offlineMolecularTransmission = permute(offlineMolecularTransmission,[1 3 2]);
offlineCombinedTransmission = permute(offlineCombinedTransmission,[1 3 2]);

int = [1:50 size(offlineMolecularTransmission,3)-50:size(offlineMolecularTransmission,3)];
eta_m=(mean(offlineCombinedTransmission(:,:,int)./offlineMolecularTransmission(:,:,int))+1).^-1;

 % etam = (mean(CalInfo.ScanData.O2OnlineComb.Transmission./CalInfo.ScanData.O2OnlineMol.Transmission)+1)^-1;
 % eta_m = (mean(onlineCombinedTransmission./onlineMolecularTransmission)+1)^-1;
 eta_c = 1-eta_m;

 etam = eta_m;
 etac = eta_c;


% %normalize

offlineMolecularTransmission = offlineMolecularTransmission./(max(offlineCombinedTransmission) .*etam./etac);
offlineCombinedTransmission  = offlineCombinedTransmission ./max(offlineCombinedTransmission);

aerosolSpectrum = zeros(size(Spectrum.nu_scan_3D_short));
aerosolSpectrum(:,:,Spectrum.offline_index(1))=1;
aerosolSpectrum = aerosolSpectrum./trapz(permute(Spectrum.nu_scan_3D_short,[3 2 1]),aerosolSpectrum,3);

Cmm = trapz(permute(Spectrum.nu_scan_3D_short,[3 2 1]),sponts6_off.*offlineMolecularTransmission,3);
Cmc = trapz(permute(Spectrum.nu_scan_3D_short,[3 2 1]),sponts6_off.*offlineCombinedTransmission,3);
Cam = trapz(permute(Spectrum.nu_scan_3D_short,[3 2 1]),aerosolSpectrum.*offlineMolecularTransmission,3);
Cac = trapz(permute(Spectrum.nu_scan_3D_short,[3 2 1]),aerosolSpectrum.*offlineCombinedTransmission,3);

Cmm = Cmm./Cac;
Cam = Cam./Cac;
Cmc = Cmc./Cac;
Cac = Cac./Cac;

Bm = 5.45e-32 .* Pr.*101325./Tr./1.380649e-23.*550^4./(lambdaCentral.^4);
Counts.mol = Data.MCS.Channel0.Data-mean(Data.MCS.Channel0.Data(470:end,:),1);
Counts.Comb = Data.MCS.Channel2.Data-mean(Data.MCS.Channel2.Data(470:end,:),1);
Gm=1;
Sm = (Gm.*Counts.mol-Cam.*Counts.Comb)./(Cmm-Cam);
Ba = (Counts.Comb./Sm-1).*Bm;


BSR = Ba./Bm+1;
%%
figure
imagesc(Time,Range.rm_raw_o2,BSR)
clim([1 8])
colormap('hot')
set(gca,'YDir','normal')
colorbar


%overlap appox
Sm = (8/3)*pi;
Tm = exp(-37.5*cumtrapz(Bm*Sm));
figure
plot(Range.rm_raw_o2,Counts.mol.*Range.rm_raw_o2.^2./Bm./(Tm.^2))
ylim([-1 5]*10^16)
xlim([0 6000])