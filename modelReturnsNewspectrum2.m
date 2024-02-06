%model of returns

clear all
%Pulse length model
h = 6.62607004e-34;
c = 2.9989e8;   % Speed of light m/s

%Input Variables
dt = .5;         %pulse duration in microseconds
nt = 40;        %number of points to calculate in each dt. THIS MUST BE AN EVEN NUMBER
nt = 2;
nt = 20;

Ts = 273;       %Surface Temperature in K      
Ps = .825;         %Surface Pressure in atm

lambda_online = 769.7958;
nu_online = 10^7./lambda_online;                    %[cm-1] Online wavenumber
WV=0;

% Range Vector
dr = c*dt*1e-6/nt;
rm = dr:dr:18000;             %Range vector in m
rkm =rm./1000;               %Range vector in km
size_r = length(rm);

% Temperature and Pressure Profiles
T = Ts-rkm.*6.5;           %Temperature profile in K
P = Ps.*(Ts./T).^-5.2199;   %Pressure profile in atm
T = T';
P = P';

%Radiosonde proflies
sondepath = 'C:\Users\Owen\OneDrive - Montana State University\Research\O2 DIAL\Data\MSU data\Radiosondes';
span_days = datetime(2022,6,22,'TimeZone','UTC');%yyyy,mm,dd)
[sonde_datetime,sondeStruc] =  COBradiosonde(sondepath,span_days);

sondeT = sondeStruc(1).T;
sondeP = sondeStruc(1).P.*0.000986923;
sondeH = sondeStruc(1).Height-sondeStruc(1).Height(1);

% % T = interp1(sondeH,sondeT,rm)';
% % P = interp1(sondeH,sondeP,rm)';

A = 1;%area of telescope
eta_O = 1;%overlap
eta_D = .6*(1/18);%detector
eta_R = 7.4516e-22 * 1.4122;%reciever
eta_R = .01;%reciever
E_pulse_on = 5 * 10^-6;%pulse energy (J)
E_pulse_off = 5 * 10^-6;%pulse energy (J)

%--Outgoint
pulse_rate = 7000;%(Hz)
avg_time = 2;%(sec)

%%

%---Backscatter
Bm = 374280*101.325*P./T./(lambda_online.^4);
% 
% Ba = Bm;
% 
% Ba = zeros(size(Bm)).*0.01;
% Ba(1:512) = 3.*Bm(1:512) ;
% %Ba(512:600) = 2.*Bm(512:600) -2.*Bm(512:600).*(1:89)'.*(1/89);
% 
% Ba(512:700) = 3.*Bm(512:700) -3.*Bm(512:700).*(1:189)'.*(1/189);
%%
BSR=[
NaN
NaN
NaN
NaN
NaN
3.33348686081786
3.36551751422900
3.24830120190958
3.51397407672946
3.24886576227332
3.28256129485638
3.20222384944337
2.99930868312942
2.86916198265844
2.64385031814741
2.47578061846368
2.42595882010490
2.36780471865549
2.44772042597672
2.37986315656486
2.40141244360314
2.52147184604419
2.65145044377840
2.63842063670806
2.70310582286966
2.72762629409563
2.67901534569066
2.55623055902631
2.53993907663394
2.49624738449374
2.69396749176289
2.90236719722746
3.13748868430628
3.10623558305595
3.28708912417195
3.02949145827879
2.40922813682599
1.76943604070225
1.41563749633709
1.19546276867337
1.10688539707404
1.11987369923108
1.09801303504226
1.13843123587360
1.11719262951284
1.10561191716656
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN];

Ba=[
    NaN
NaN
NaN
NaN
NaN
6.32907992830923e-07
6.37710590617329e-07
6.02431327959596e-07
6.69518780572237e-07
5.95261370654154e-07
6.00484865622041e-07
5.75797461542331e-07
5.19529228403471e-07
4.82716417259260e-07
4.21905813186831e-07
3.76423065828719e-07
3.61456164568666e-07
3.44556130881516e-07
3.62410512957459e-07
3.43261612072238e-07
3.46434565652825e-07
3.73747293067504e-07
4.03117256846504e-07
3.97407212971547e-07
4.10477358414255e-07
4.13739821262410e-07
3.99534967784686e-07
3.67950426831804e-07
3.61765163599222e-07
3.49242422286174e-07
3.92845292441800e-07
4.38324953152219e-07
4.89309294805676e-07
4.79023446621442e-07
5.16767802649679e-07
4.55569274199370e-07
3.14264704622478e-07
1.70461440408553e-07
9.14744603992016e-08
4.27339292735071e-08
2.32136212586037e-08
2.58616491946922e-08
2.10046738413474e-08
2.94684951760562e-08
2.47803728997268e-08
2.21817605779635e-08
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN];
Bm=[
2.79528095465145e-07
2.77856427388800e-07
2.76190633539672e-07
2.74530707835691e-07
2.72876644185929e-07
2.71228436490573e-07
2.69586078640885e-07
2.67949564519169e-07
2.66318887998735e-07
2.64694042943861e-07
2.63075023209760e-07
2.61461822642539e-07
2.59854435079168e-07
2.58252854347438e-07
2.56657074265928e-07
2.55067088643964e-07
2.53482891281583e-07
2.51904475969496e-07
2.50331836489049e-07
2.48764966612182e-07
2.47203860101395e-07
2.45648510709708e-07
2.44098912180616e-07
2.42555058248059e-07
2.41016942636376e-07
2.39484559060264e-07
2.37957901224742e-07
2.36436962825110e-07
2.34921737546905e-07
2.33412219065860e-07
2.31908401047868e-07
2.30410277148934e-07
2.28917841015135e-07
2.27431086282582e-07
2.25950006577369e-07
2.24474595515538e-07
2.23004846703032e-07
2.21540753735652e-07
2.20082310199014e-07
2.18629509668504e-07
2.17182345709234e-07
2.15740811875998e-07
2.14304901713227e-07
2.12874608754940e-07
2.11449926524705e-07
2.10030848535588e-07
2.08617368290107e-07
2.07209479280188e-07
2.05807174987115e-07
2.04410448881485e-07
2.03019294423160e-07
2.01633705061216e-07
2.00253674233901e-07
1.98879195368578e-07
1.97510261881683e-07
1.96146867178671e-07
1.94789004653969e-07
1.93436667690925e-07
1.92089849661757e-07
1.90748543927500e-07
1.89412743837961e-07
1.88082442731659e-07
1.86757633935780e-07
1.85438310766120e-07
1.84124466527035e-07
1.82816094511386e-07
1.81513188000485e-07
1.80215740264043e-07
1.78923744560112e-07
1.77637194135034e-07
1.76356082223384e-07
1.75080402047912e-07
1.73810146819491e-07
1.72545309737057e-07
1.71285883987552e-07
1.70031862745871e-07
1.68783239174797e-07
1.67540006424946e-07
1.66302157634712e-07
1.65069685930198e-07];

hsrlRange=[
    22.4221717500000
97.3702862500000
172.318400750000
247.266515250000
322.214629750000
397.162744250000
472.110858750000
547.058973250000
622.007087750000
696.955202250000
771.903316750000
846.851431250000
921.799545750000
996.747660250000
1071.69577475000
1146.64388925000
1221.59200375000
1296.54011825000
1371.48823275000
1446.43634725000
1521.38446175000
1596.33257625000
1671.28069075000
1746.22880525000
1821.17691975000
1896.12503425000
1971.07314875000
2046.02126325000
2120.96937775000
2195.91749225000
2270.86560675000
2345.81372125000
2420.76183575000
2495.70995025000
2570.65806475000
2645.60617925000
2720.55429375000
2795.50240825000
2870.45052275000
2945.39863725000
3020.34675175000
3095.29486625000
3170.24298075000
3245.19109525000
3320.13920975000
3395.08732425000
3470.03543875000
3544.98355325000
3619.93166775000
3694.87978225000
3769.82789675000
3844.77601125000
3919.72412575000
3994.67224025000
4069.62035475000
4144.56846925000
4219.51658375000
4294.46469825000
4369.41281275000
4444.36092725000
4519.30904175000
4594.25715625000
4669.20527075000
4744.15338525000
4819.10149975000
4894.04961425000
4968.99772875000
5043.94584325000
5118.89395775000
5193.84207225000
5268.79018675000
5343.73830125000
5418.68641575000
5493.63453025000
5568.58264475000
5643.53075925000
5718.47887375000
5793.42698825000
5868.37510275000
5943.32321725000];
%%

Ba = fillmissing(interp1(hsrlRange,Ba',rm),'nearest');
Bm = fillmissing(interp1(hsrlRange,Bm',rm),'nearest');

Ba = smoothdata(Ba,'movmean',20);
Bm = smoothdata(Bm,'movmean',20);

Ba = Ba';
Bm = Bm';

%Ba=zeros(size(Bm));
Ba = Bm*1;

 Ba(1:268)=Bm(1:268)*2;
 Ba(269:(269+149)) = Bm(269:(269+149)) .* linspace(2,1,150)';

%--transmission
Sa = 60;
Ta = exp(-dr*cumtrapz(Ba*Sa));
Sm = (8/3)*pi;
Tm = exp(-dr*cumtrapz(Bm*Sm));

%%
%--spectrum
Spectrum.lambda_online = 769.7958;
Spectrum.lambda_offline = 770.1085;

Spectrum.nu_online = 10^7./Spectrum.lambda_online;                    %[cm-1] Online wavenumber
Spectrum.nu_offline = 10^7./Spectrum.lambda_offline;                  %[cm-1] Offline wavenumber

nuMin = Spectrum.nu_online-0.334*5;                                 %[cm-1] Scan lower bound
nuMax = Spectrum.nu_online+0.334*5;                                 %[cm-1] Scan upper bound
nuMin = Spectrum.nu_online-0.334;                                 %[cm-1] Scan lower bound
nuMax = Spectrum.nu_online+0.334;  
Spectrum.nuBin = 0.00222;                                    %[cm-1] Scan increment
%Spectrum.nuBin = 0.00222/3;  
%Spectrum.nuBin = 4.6667e-4;
nu_scan = (nuMin:Spectrum.nuBin:nuMax);                      %[cm-1](1 x nu) Scan vector

nuMin_off = Spectrum.nu_offline-0.334*5;                                 %[cm-1] Scan lower bound
nuMax_off = Spectrum.nu_offline+0.334*5;                                 %[cm-1] Scan upper bound
nuMin_off = Spectrum.nu_offline-0.334;                                 %[cm-1] Scan lower bound
nuMax_off = Spectrum.nu_offline+0.334; 
nu_scan_off = (nuMin_off:Spectrum.nuBin:nuMax_off);

Spectrum.nu_scan_3D_short = permute(nu_scan, [3 1 2]);       %[cm-1] putting scan in third dimension
Spectrum.nu_scan_3D_short_off = permute(nu_scan_off, [3 1 2]);       %[cm-1] putting scan in third dimension
Spectrum.lambda_scan_3D_short_off = 10^7./Spectrum.nu_scan_3D_short_off;

[~,online_index] = min(abs(Spectrum.nu_online-Spectrum.nu_scan_3D_short));
[~,offline_index] = min(abs(Spectrum.nu_offline-Spectrum.nu_scan_3D_short_off));

%--O2absorption
%%%o2absorption = absorption_O2_770_model_wavenumber(T,P,Spectrum.nu_scan_3D_short,WV);
o2absorption = absorption_O2_770_PCA(T,P,Spectrum.nu_scan_3D_short,WV);
%o2absorption_off = absorption_O2_770_model_wavenumber(T,P,Spectrum.nu_scan_3D_short_off ,WV);
%o2absorption_off = absorption_O2_770_PCA_off(T,P,Spectrum.nu_scan_3D_short_off,WV);
o2absorption_off = absorption_O2_770_model(T,P,Spectrum.nu_offline,WV);
absorption = absorption_O2_770_model(T,P,nu_online,WV);
TO2 = exp(-dr*cumtrapz(absorption));

TO2nu = exp(-dr*cumtrapz(o2absorption));
TO2nu_off = exp(-dr*cumtrapz(o2absorption_off));

%%

%--molecularbackscatter

    cB = 1.2;%Brullouion correction to doppler gaussian half width
    cB = -0.01.*((1500+rm')/1000) + 1.2;
    m_air = 28.97/1000./6.02214e23;
    kb = 1.3806e-23;
    c=3e8;

    c_doppler_O2 = m_air*c^2./(8*(Spectrum.nu_online(1)*100).^2*kb);                   %[m^2 K] Doppler coefficient
    doppler_O2_un_ret = ((c_doppler_O2./T/pi).^0.5).*exp(-c_doppler_O2.*(Spectrum.nu_online(1)*100-Spectrum.nu_scan_3D_short*100).^2./T./cB.^2); %[m] Doppler broadended lineshape         

    norm_O2_ret = trapz(doppler_O2_un_ret,3).*Spectrum.nuBin*100;                   %[none] Lineshape integral
    doppler_O2_ret = doppler_O2_un_ret./norm_O2_ret;                       %[m] Normalized doppler lineshape

    c_doppler_O2 = m_air*c^2./(8*(Spectrum.nu_offline(1)*100).^2*kb);                   %[m^2 K] Doppler coefficient
    doppler_O2_un_ret = ((c_doppler_O2./T/pi).^0.5).*exp(-c_doppler_O2.*(Spectrum.nu_offline(1)*100-Spectrum.nu_scan_3D_short_off*100).^2./T./cB.^2); %[m] Doppler broadended lineshape         

    norm_O2_ret = trapz(doppler_O2_un_ret,3).*Spectrum.nuBin*100;                   %[none] Lineshape integral
    doppler_O2_ret_off = doppler_O2_un_ret./norm_O2_ret;                       %[m] Normalized doppler lineshape

gaerosol = zeros(size(Spectrum.nu_scan_3D_short));
gaerosol(online_index) = 1./(Spectrum.nuBin*100);

gaerosol_off = zeros(size(Spectrum.nu_scan_3D_short_off));
gaerosol_off(offline_index) = 1./(Spectrum.nuBin*100);

[Spectrum] = PCAconstrunctionRB2(Spectrum);
[sponts6] = RB_O2_770_PCA(T,P,Spectrum.nu_scan_3D_short,Spectrum);
[sponts6_off] = RB_O2_770_PCA(T,P,Spectrum.nu_scan_3D_short_off,Spectrum);
%[sponts6_off] = RB_O2_770_PCA_offline(T,P,Spectrum.nu_scan_3D_short_off,Spectrum);
doppler_O2_ret = sponts6;
doppler_O2_ret_off = sponts6_off;


% gaerosol = zeros(size(Spectrum.nu_scan_3D_short));
% gaerosol(online_index:online_index+2) = 1./(Spectrum.nuBin*100)/3;
% 
% gaerosol_off = zeros(size(Spectrum.nu_scan_3D_short_off));
% gaerosol_off(offline_index:offline_index+2) = 1./(Spectrum.nuBin*100)/3;

transmitPhotons_on = E_pulse_on./(h*c./(lambda_online*10^-9)) * (pulse_rate/2) * avg_time*60 /10;
transmitPhotons_off = E_pulse_off./(h*c./(Spectrum.lambda_offline*10^-9)) * (pulse_rate/2) * avg_time*60 /10;
%%
%%
%spectral purity
% % spectralpurity=1;
% % aseWidth = 10e9;%hz
% % 
% % aseWidth = 100e6;%hz
% % %aseWidth = 0;
% % ase=normpdf(Spectrum.nu_scan_3D_short,Spectrum.nu_online,aseWidth./c./100);
% % ase = ase./trapz(ase)./(Spectrum.nuBin*100);
% % ase_off=normpdf(Spectrum.nu_scan_3D_short_off,Spectrum.nu_offline,aseWidth./c./100);
% % ase_off = ase_off./trapz(ase_off)./(Spectrum.nuBin*100);
% % 
% % laserWidth = gaerosol.*spectralpurity + ase.*(1-spectralpurity);
% % laserWidth_off = gaerosol_off.*spectralpurity + ase_off.*(1-spectralpurity);

spectralpurity =.999;
%spectralpurity =1;
gaerosol(online_index) = gaerosol(online_index).*spectralpurity;
gaerosol(1:online_index-1) = gaerosol(online_index).*(1-spectralpurity)./(length(gaerosol)-1);
gaerosol(online_index+1:end) = gaerosol(online_index).*(1-spectralpurity)./(length(gaerosol)-1);

gaerosol_off(offline_index) = gaerosol_off(offline_index).*spectralpurity;
gaerosol_off(1:offline_index-1) = gaerosol_off(offline_index).*(1-spectralpurity)./(length(gaerosol)-1);
gaerosol_off(offline_index+1:end) = gaerosol_off(offline_index).*(1-spectralpurity)./(length(gaerosol)-1);

laserWidth = gaerosol;

laserWidth_off = gaerosol_off;

aseWidth = ones(size(laserWidth))./(Spectrum.nuBin*100)./size(laserWidth,3 );

%%
%Overlap

nsPerBin = double(250);                             %[ns] Nanoseconds per bin
NBins = double(560);                                %[none] Number of range bins
rangeBin = (c * nsPerBin(1)*10^-9)/2;               %[m] Create range bin size from speed of light and bin time over 2
rangeMin = -150; %[m]
rangeMin = -300; %[m]
rangeMin = -rangeBin; %[m]
rangeMin = -(c * (1*10^-6))/2;
rangeMin = 0;
rm_raw_o2 = rangeMin:rangeBin:NBins(1)*rangeBin+rangeMin-rangeBin;    %[m] Create range vector
rm_raw_o2 = rm_raw_o2(:);                           %[m] Convert range vector to column vector
r_max = 6000;   %[m] Max range 
%r_max = 10000;   %[m] Max range 
%%
rm_over = rm_raw_o2(rm_raw_o2<=r_max & rm_raw_o2>0);     %[m] Shorten range vector to max range
%load('D:\OneDrive - Montana State University\Research\ARM data analysis\overlap7_21_20.mat','overlap')
% load('C:\Users\Owen\OneDrive - Montana State University\Research\ARM data analysis\overlap7_21_20.mat','overlap')

%eta_O = interp1(rm_over,overlap,rm);
%eta_O = fillmissing(eta_O,'nearest');

  load('OverlapSim_5_31_23.mat','R','OVF','OVF_near')
  eta_O = interp1(R,OVF,rm);
  eta_O = fillmissing(eta_O,'nearest');

eta_O_near = interp1(R,OVF_near,rm);
 eta_O_near = fillmissing(eta_O_near,'nearest');


   load('OverlapSimASE_7_12_23.mat','R','OVF','OVF_near')
   %load('OverlapSimASE_7_35_23.mat','R','OVF','OVF_near')
  eta_OASE = interp1(R,OVF,rm);
  eta_OASE = fillmissing(eta_OASE,'nearest');

eta_O_nearASE = interp1(R,OVF_near,rm);
 eta_O_nearASE = fillmissing(eta_O_nearASE,'nearest');

 %%

Fon  = eta_O'*.98.*Tm.^2.*Ta.^2.*TO2.^2.*(Bm+Ba)./rm'.^2;
Foff  = eta_O'*.98.*Tm.^2.*Ta.^2.*(Bm+Ba)./rm'.^2;

Fon  = eta_O'*.98.*Tm.^2.*Ta.^2.*TO2.^2.*(Bm+Ba)./rm'.^2 + eta_O_near'*.02.*Tm.^2.*Ta.^2.*TO2.^2.*(Bm+Ba)./rm'.^2;
Foff  = eta_O'*.98.*Tm.^2.*Ta.^2.*(Bm+Ba)./rm'.^2 + eta_O_near'*.02.*Tm.^2.*Ta.^2.*(Bm+Ba)./rm'.^2;



FonSpectrumUp = trapz(Tm.*Ta.*TO2nu.*laserWidth,3).*Spectrum.nuBin*100;
FonSpectrumDown = trapz(Tm.*Ta.*TO2nu.*(Bm.*(convn(doppler_O2_ret,laserWidth,'same').*Spectrum.nuBin.*100) ...
    + Ba.*laserWidth),3).*Spectrum.nuBin*100;

FoffSpectrumUp = trapz(Tm.*Ta.*TO2nu_off.*laserWidth_off,3).*Spectrum.nuBin*100;
FoffSpectrumDown = trapz(Tm.*Ta.*TO2nu_off.*(Bm.*(convn(doppler_O2_ret_off,laserWidth_off,'same').*Spectrum.nuBin.*100) ...
    + Ba.*laserWidth_off),3).*Spectrum.nuBin*100;

FonSpectrumUpASE = trapz(Tm.*Ta.*TO2nu.*aseWidth,3).*Spectrum.nuBin*100;
FonSpectrumDownASE = trapz(Tm.*Ta.*TO2nu.*(Bm.*(convn(doppler_O2_ret,aseWidth,'same').*Spectrum.nuBin.*100) ...
    + Ba.*aseWidth),3).*Spectrum.nuBin*100;


FoffSpectrumUpASE = trapz(Tm.*Ta.*TO2nu_off.*aseWidth,3).*Spectrum.nuBin*100;
FoffSpectrumDownASE = trapz(Tm.*Ta.*TO2nu_off.*(Bm.*(convn(doppler_O2_ret_off,aseWidth,'same').*Spectrum.nuBin.*100) ...
    + Ba.*aseWidth),3).*Spectrum.nuBin*100;

%%
N_on  = transmitPhotons_on  .* eta_R .* eta_D .* A .* dr .* Fon *2.5;
N_off = transmitPhotons_on .* eta_R .* eta_D .* A .* dr .* Foff *2.5;

load('Afterpulsing_correction_09092022.mat','Correction_Nc_on','Correction_Nc_off')
correctionRange = 0:(250e-9*3e8/2):(length(Correction_Nc_on)-1)*(250e-9*3e8/2);
Correction_Nc_on = interp1(correctionRange,Correction_Nc_on,rm)';
Correction_Nc_off = interp1(correctionRange,Correction_Nc_off,rm)';

%Correction_Nc_on = [Correction_Nc_on'; mean(Correction_Nc_on(:,end-20:end)).*ones(length(N_on)-length(Correction_Nc_on),1)];
%Correction_Nc_off = [Correction_Nc_off'; mean(Correction_Nc_off(:,end-20:end)).*ones(length(N_off)-length(Correction_Nc_off),1)];
N_onAfterpulse = N_on + movmean((Correction_Nc_on),5)-mean(Correction_Nc_on(end-20*4,end));
N_offAfterpulse = N_off + movmean((Correction_Nc_off),5)-mean(Correction_Nc_off(end-20*4,end));

N_onSpectrum = eta_O'.*.98.*FonSpectrumUp.*FonSpectrumDown.*eta_R.*eta_D.*A.*dr.*transmitPhotons_on./rm'.^2 + eta_O_near'.*.02.*FonSpectrumUp.*FonSpectrumDown.*eta_R.*eta_D.*A.*dr.*transmitPhotons_on./rm'.^2;
N_offSpectrum = eta_O'.*.98.*FoffSpectrumUp.*FoffSpectrumDown.*eta_R.*eta_D.*A.*dr.*transmitPhotons_off./rm'.^2 + eta_O_near'.*.02.*FoffSpectrumUp.*FoffSpectrumDown.*eta_R.*eta_D.*A.*dr.*transmitPhotons_off./rm'.^2;


N_onSpectrumASE = eta_OASE'.*.98.*FonSpectrumUpASE.*FonSpectrumDownASE.*eta_R.*eta_D.*A.*dr.*transmitPhotons_on*.01./rm'.^2 + eta_O_nearASE'.*.02.*FonSpectrumUpASE.*FonSpectrumDownASE.*eta_R.*eta_D.*A.*dr.*transmitPhotons_on*.01./rm'.^2;
N_offSpectrumASE = eta_OASE'.*.98.*FoffSpectrumUpASE.*FoffSpectrumDownASE.*eta_R.*eta_D.*A.*dr.*transmitPhotons_off*.01./rm'.^2 + eta_O_nearASE'.*.02.*FoffSpectrumUpASE.*FoffSpectrumDownASE.*eta_R.*eta_D.*A.*dr.*transmitPhotons_off*.01./rm'.^2;

%%

BSR = (Ba+Bm)./Bm;
smoothingLength = 300;%pulse length
smoothingLength = 150;%pulse length
smoothPoints = round(smoothingLength/dr);
%smoothPoints=1;
smoothVector = ones(smoothPoints,1)/smoothPoints;


%%% Summing over pulse length
for iii = 1:length(rm)
    if iii <=smoothPoints
        
        N_onPulse(iii,:) = mean(N_on(1:iii,:),1);
        N_offPulse(iii,:) = mean(N_off(1:iii,:),1);
    
        N_onSpectrumPulse(iii,:) = mean(N_onSpectrum(1:iii,:),1);
        N_offSpectrumPulse(iii,:) = mean(N_offSpectrum(1:iii,:),1);

        BSRshift(iii,:) = mean(BSR(1:iii,:),1);


    else
        N_onPulse(iii,:) = sum(N_on(iii-smoothPoints:iii,:),1)./smoothPoints;
        N_offPulse(iii,:) = sum(N_off(iii-smoothPoints:iii,:),1)./smoothPoints;
    
        N_onSpectrumPulse(iii,:) = sum(N_onSpectrum(iii-smoothPoints:iii,:),1)./smoothPoints;
        N_offSpectrumPulse(iii,:) = sum(N_offSpectrum(iii-smoothPoints:iii,:),1)./smoothPoints;

        BSRshift(iii,:) = sum(BSR(iii-smoothPoints:iii,:),1)./smoothPoints;
        
    end
end

N_onPulseAfterpulse = N_onPulse- movmean((Correction_Nc_on),5)-mean(Correction_Nc_on(end-20*4,end));
N_offPulseAfterpulse = N_offPulse- movmean((Correction_Nc_off),5)-mean(Correction_Nc_off(end-20*4,end));


% BSRshift = [ones(smoothPoints,1).*BSR(1); BSR(1:(end-smoothPoints))];

%%
[a_0] = alpha_0(N_on,N_off,dr);
[a_0Pulse] = alpha_0(N_onPulse,N_offPulse,dr);
[a_0Spectrum] = alpha_0(N_onSpectrum,N_offSpectrum,dr);
[a_0SpectrumPulse] = alpha_0(N_onSpectrumPulse,N_offSpectrumPulse,dr);
[a_0Afterpulse] = alpha_0(N_onAfterpulse,N_offAfterpulse,dr);
[a_0PulseAfterpulse] = alpha_0(N_onPulseAfterpulse,N_offPulseAfterpulse,dr);


alpha_0 = a_0+o2absorption_off;

figure()
plot(a_0,rm)
hold on
plot(a_0Pulse,rm)
plot(a_0Spectrum,rm)
plot(a_0SpectrumPulse,rm)
plot(a_0Afterpulse,rm)
plot(a_0PulseAfterpulse,rm)
plot(o2absorption(:,:,online_index),rm,'--')

legend('a_0','a_0Pulse','a_0Spectrum','a_0SpectrumPulse','a_0 afterpulse','absorption')

figure()
plot(a_0-o2absorption(:,:,online_index),rm)
hold on
plot(a_0Pulse-o2absorption(:,:,online_index),rm)
plot(a_0Spectrum-o2absorption(:,:,online_index),rm)
plot(a_0SpectrumPulse-o2absorption(:,:,online_index),rm)
plot(a_0Afterpulse-o2absorption(:,:,online_index),rm)
plot(a_0PulseAfterpulse-o2absorption(:,:,online_index),rm)

legend('a_0','a_0Pulse','a_0Spectrum','a_0SpectrumPulse','a_0 afterpulse','absorption')

%%

%--correction
altitude = 1.5719;
Model.T = T;
Model.P = P;
Model.WV = 0;
BSR = (Ba+Bm)./Bm;

Options.oversample=1;
Range.rm = rm';
Range.i_range = length(rm);
Range.rangeBin = rm(2)-rm(1);
Time.ts = 1;
Time.i_time = 1;
Spectrum.i_scan_3D_short = length(Spectrum.nu_scan_3D_short);
Spectrum.online_index = online_index;
Spectrum.offline_index = offline_index;
[alpha_final,alpha_1_raw,alpha_2_raw,Spectrum] = pertAbsorption(a_0, 1, Model, Range, Time, Spectrum, BSR, 0,0, Options, 0,false);
Alpha_total = a_0+alpha_1_raw+alpha_2_raw;

[alpha_final,alpha_1_raw,alpha_2_raw,Spectrum] = pertAbsorption(a_0Pulse, 1, Model, Range, Time, Spectrum, BSR, 0,0, Options, 0,false);
Alpha_totalPulse = a_0Pulse+alpha_1_raw+alpha_2_raw;

[alpha_final,alpha_1_raw,alpha_2_raw,Spectrum] = pertAbsorption(a_0Spectrum, 1, Model, Range, Time, Spectrum, BSR, 0,0, Options, 0,false);
Alpha_totalSpectrum = a_0Spectrum+alpha_1_raw+alpha_2_raw;

[alpha_final,alpha_1_raw,alpha_2_raw,Spectrum] = pertAbsorption(a_0SpectrumPulse, 1, Model, Range, Time, Spectrum, BSR, 0,0, Options, 0,false);
Alpha_totalSpectrumPulse = a_0SpectrumPulse+alpha_1_raw+alpha_2_raw;

[alpha_final,alpha_1_raw,alpha_2_raw,Spectrum] = pertAbsorption(a_0SpectrumPulse, 1, Model, Range, Time, Spectrum, BSRshift, 0,0, Options, 0,false);
Alpha_totalSpectrumPulseBSRShift = a_0SpectrumPulse+alpha_1_raw+alpha_2_raw;

figure()
plot(Alpha_total,rm)
hold on
plot(Alpha_totalPulse,rm)
plot(Alpha_totalSpectrum,rm)
plot(Alpha_totalSpectrumPulse,rm)
plot(Alpha_totalSpectrumPulseBSRShift,rm)
plot(o2absorption(:,:,online_index),rm,'--')
legend('a_t','a_tPulse','a_tSpectrum','a_tSpectrumPulse','absorption')

%%

[T_final,Lapse,Ts_fit,P_final,mean_lapse_rate,exclusion,Titer] =  temperatureRetrieval(T,1,rm',0,WV,Spectrum.nu_online,a_0,0,logical(zeros(size(rm'))),T(1),P(1),-6.5/1000);
[T_finalPulse,Lapse,Ts_fit,P_final,mean_lapse_rate,exclusion,Titer] =  temperatureRetrieval(T,1,rm',0,WV,Spectrum.nu_online,a_0Pulse,0,logical(zeros(size(rm'))),T(1),P(1),-6.5/1000);
[T_finalSpectrum,Lapse,Ts_fit,P_final,mean_lapse_rate,exclusion,Titer] =  temperatureRetrieval(T,1,rm',0,WV,Spectrum.nu_online,Alpha_totalSpectrum,0,logical(zeros(size(rm'))),T(1),P(1),-6.5/1000);
[T_finalSpectrumPulse,Lapse,Ts_fit,P_final,mean_lapse_rate,exclusion,Titer] =  temperatureRetrieval(T,1,rm',0,WV,Spectrum.nu_online,Alpha_totalSpectrumPulse,0,logical(zeros(size(rm'))),T(1),P(1),-6.5/1000);
[T_finalSpectrumPulseBSRShift,Lapse,Ts_fit,P_final,mean_lapse_rate,exclusion,Titer] =  temperatureRetrieval(T,1,rm',0,WV,Spectrum.nu_online,Alpha_totalSpectrumPulseBSRShift,0,logical(zeros(size(rm'))),T(1),P(1),-6.5/1000);

figure
plot(T_final,rm)
hold on
plot(T_finalPulse,rm)
plot(T_finalSpectrum,rm)
plot(T_finalSpectrumPulse,rm)
plot(T_finalSpectrumPulseBSRShift,rm)
plot(T,rm,'--')
legend('Tfinal','TfinalPulse','TfinalSpectrum','TfinalSpectrumPulse','shift','T')

figure
plot(T-T_final,rm)
hold on
plot(T-T_finalPulse,rm)
plot(T-T_finalSpectrum,rm)
plot(T-T_finalSpectrumPulse,rm)
plot(T-T_finalSpectrumPulseBSRShift,rm,'--')
plot(T(1:end-smoothPoints)-T_finalSpectrumPulseBSRShift((smoothPoints+1):end),rm(1:end-smoothPoints),'--')
legend('Tfinal','TfinalPulse','TfinalSpectrum','TfinalSpectrumPulse')


