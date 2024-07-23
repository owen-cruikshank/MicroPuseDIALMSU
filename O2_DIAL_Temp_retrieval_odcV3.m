% O2_DIAL_Temp_retrieval_odcV3.m
% Analysis program for O2 DIAL instrument from 
clear all

%=Date Range
date_start = datetime(2022,4,21,'TimeZone','UTC');%yyyy,mm,dd
date_end = datetime(2022,7,16,'TimeZone','UTC');%yyyy,mm,dd

date_start = datetime(2022,7,27,'TimeZone','UTC');%yyyy,mm,dd
date_end = datetime(2022,09,12,'TimeZone','UTC');%yyyy,mm,dd

date_start = datetime(2024,07,12,'TimeZone','UTC');%yyyy,mm,dd
date_start = datetime(2024,07,19,'TimeZone','UTC');%yyyy,mm,dd
date_end = datetime(2024,07,19,'TimeZone','UTC');%yyyy,mm,dd

span_days = date_start:date_end;

%=Time and range averaging
Options.intTime = 10;  %[min] Integration time
Options.intRange = 2; %[bins] Integration range

Options.t_avg = 1;     %[bins] Time smoothing bins
Options.oversample = 1; %[bins] Range smoothing bins

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

%======================
%==== Load data
%==== photon counts, radiosondes, weather station data, HSRL data, etc.
%======================

Options.MPDname = '00';
%Options.MPDname = '03';
%Options.MPDname = '05';

Options.DataPath = 'C:\Users\Owen\OneDrive - Montana State University\Research\O2 DIAL\Data';
if strcmp(Options.MPDname,'00')
    %==Load data from MSU instument==
    [Range,Time,Counts,Sonde,Model,Spectrum,HSRL,Data] = loadMSUdata(span_days,Options,Constants);

% ==Load data from SGP ARM instument==
%date_begin = datetime(2019,4,17); date_end   = datetime(2019,4,22);
%span_days = date_begin:date_end;        % Days in set [datetime]
%[Range,Time,Counts,Sonde,Model,Spectrum,BSR,Data] = loadSGPdata(span_days,Options,Constants);

%==Load data from Boulder instument==
elseif strcmp(Options.MPDname,'03') || strcmp(Options.MPDname,'01') || strcmp(Options.MPDname,'05')
    [Range,Time,Counts,Sonde,Model,Spectrum,HSRL,Data] =loadBoulderdata3(span_days,Options,Constants);
end
%%
%======================
%= Cloud and SNR mask =
%======================
disp('Calculating masks')
 
cloud_p_point = 200; %(hours) Time to plot mask data

%Threshold for initial SNR mask
SNR_threshold = 2300;

%Standard deviation Threshold for Cloud mask
SD_threshold = 5;

BGmult =1; % Multiplier for background for SNR calculation
%Low altitude masking
lowAlt = 330;
[SNRm , cloud_SDm_above, cloud_SDm,~] = mask_O2_BSR(cloud_p_point,SNR_threshold,SD_threshold,Options.oversample,Options.t_avg,Counts,BGmult,Time,Range,HSRL.BSR,lowAlt);

%Create initial mask logicals above clouds
cloud_SDm_above(cloud_SDm_above==-1)=0;
cloud_SDm_above = ~(logical(cloud_SDm_above));
SNRm = ~(logical(SNRm));
cloud_SDm_above = cloud_SDm_above | SNRm;
cloud_SDm = logical(cloud_SDm);


%%

%======================
%==== Poisson Thin
%======================
%Number of poisson thinning/bootstrapping iterations
iter = 30;
Counts = poissonThin2(Counts,cloud_SDm_above,iter);

%%
%======================
%==== Bootsrapping iteration
%======================

for jjjj = 1:iter
    display(['Bootstrapping iteration ', num2str(jjjj)])
    if jjjj >=2 %=== Set model to retrieval data ===
        % Model.T = real(fillmissing(Temperature.T_final_testFull(:,:,jjjj-1),'linear'));
        %Model.P = real(fillmissing(Temperature.Patm_finalFull,'linear'));
        Model.T = fillmissing(Temperature.L_fit_sm_test(:,:,end).*Range.rm+Temperature.Ts_fit(:,:,end),'linear');
        Model.T(Model.T<0)=0.001;
        Model.Ts =Temperature.Ts_fit(:,:,end);
        Model.P = real(fillmissing(Temperature.Patm_finalFull,'linear'));

        %Testing random lapse
        % % LapseRand = normrnd(-6.5,4);
        % % LapseRand = rand.*-5.8 -4;
        % % LapseRand = LapseRand.*ones(size(Model.Ts));
        % % Model.T = fillmissing(LapseRand.*Range.rkm+Model.Ts,'linear');

        g0 = 9.80665;               %[m/s/s]Gravitational acceleration
        M_air = 0.0289644;          %[kg/mol] Molar mass of air
        q_O2 = .2095;               %[unitless] O2 atmospheric mixing ratio 
        kB = 1.38065E-23;           %[J/K] Boltzman's constant 
        N_A = 6.02214076e23;        %[1/mol] Avagadro's number
        R = kB * N_A;               %[J/K/mol] universal gas constant
        gamma = g0 * M_air / R;     %[K/m]gravity molar mass of air and gas constant
        %Pg = fillmissing(Model.Ps.*(Temperature.Ts_fit(:,:,end)./(Temperature.Ts_fit(:,:,end)+Temperature.L_fit_sm_test(:,:,end).*Range.rm)).^(gamma./Temperature.L_fit_sm_test(:,:,end)),'linear');
      
        %Pg = fillmissing(Model.Ps.*(Model.Ts./(Model.Ts+LapseRand.*Range.rkm)).^(gamma./LapseRand/1000),'linear');
        %%%%Model.P = Pg;

        Model.absorption = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_online,Model.WV,Constants);

        Model.absorption_off = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_offline,Model.WV,Constants); %[m-1] Funcrtion to calculate theoretical absorption
        Model.transmission = exp(-cumtrapz(Range.rm,Model.absorption));
    end


    %=== Calculate HSRL ===
    if strcmp(Options.MPDname,'00')
        Atmosphere.Pressure = Model.P./0.009869233; 
        Atmosphere.Temperature = Model.T;
        Counts.Nc_on = Counts.o2on;%.*Counts.NBins;
        Counts.Nc_off = Counts.o2off;%.*Counts.NBins;
        Counts.Nm_on = Counts.o2on_mol;%.*Counts.NBins;
        Counts.Nm_off = Counts.o2off_mol;%.*Counts.NBins;
        Counts.sigma_Nm_off = sqrt(Counts.Nc_on);
        Counts.sigma_Nm_on = sqrt(Counts.Nc_off);
        Counts.sigma_Nc_off = sqrt(Counts.Nm_on);
        Counts.sigma_Nc_on = sqrt(Counts.Nm_off);
        Options.t_step = 1;
        %[HSRL] = HSRL_retrieval_20220909(Counts,Atmosphere,Options,Spectrum);
        [HSRL] = backscatterRetrievalMPD(Counts, Model, Spectrum, Options);
        %[HSRL] = HSRL_retrieval_20230115(Counts,Atmosphere,Options);

        Counts.Nc_on = Counts.o2on_noise;%.*Counts.NBins;
        Counts.Nc_off = Counts.o2off_noise;%.*Counts.NBins;
        Counts.Nm_on = Counts.o2on_noise_mol;%.*Counts.NBins;
        Counts.Nm_off = Counts.o2off_noise_mol;%.*Counts.NBins;
        [HSRLfull] = backscatterRetrievalMPD(Counts, Model, Spectrum, Options);

        %Smoothing HSRL
        k = ones(1,1)./(1*1);
        HSRLfull.BSR = nanconv(HSRLfull.BSR,k,'edge','nanout');
    
    elseif strcmp(Options.MPDname,'03') || strcmp(Options.MPDname,'05')
        Counts.Nc_on = Counts.o2on;%.*Counts.NBins;
        Counts.Nc_off = Counts.o2off;%.*Counts.NBins;
        Counts.Nm_on = Counts.o2on_mol;%.*Counts.NBins;
        Counts.Nm_off = Counts.o2off_mol;%.*Counts.NBins;
        [HSRL] = backscatterRetrievalMPD(Counts, Model, Spectrum,Options);
        % Counts.Nc_on = Counts.o2on_bgsub;%.*Counts.NBins;
        % Counts.Nc_off = Counts.o2off_bgsub;%.*Counts.NBins;
        % Counts.Nm_on = Counts.o2on_bgsub_mol;%.*Counts.NBins;
        % Counts.Nm_off = Counts.o2off_bgsub_mol;%.*Counts.NBins;

        Counts.Nc_on = Counts.o2on_noise;%.*Counts.NBins;
        Counts.Nc_off = Counts.o2off_noise;%.*Counts.NBins;
        Counts.Nm_on = Counts.o2on_noise_mol;%.*Counts.NBins;
        Counts.Nm_off = Counts.o2off_noise_mol;%.*Counts.NBins;
        [HSRLfull] = backscatterRetrievalMPD(Counts, Model, Spectrum, Options);

        %Smoothing HSRL
        k = ones(1,1)./(1*1);
        HSRL.BSR = nanconv(HSRL.BSR,k,'edge','nanout');
    end
    
    %Calulate HSRL at 828nm
    HSRL.Bm828 = HSRL.Bm *(770/828)^4;
    HSRL.Ba828 = HSRL.Ba*(770/828);
    HSRL.BSR828 = HSRL.Ba828./HSRL.Bm828+1;
  
%%
    %Calculate Bootsrapping HSRL f and g
    if strcmp(Options.MPDname,'00')

        Atmosphere.Pressure = Model.P./0.009869233; 
        Atmosphere.Temperature = Model.T;
        Counts.Nc_on = Counts.fon(:,:,jjjj);
        Counts.Nc_off = Counts.foff(:,:,jjjj);
        Counts.Nm_on = Counts.fon_mol(:,:,jjjj);
        Counts.Nm_off = Counts.foff_mol(:,:,jjjj);
        Counts.sigma_Nm_off = 0;
        Counts.sigma_Nm_on = 0;
        Counts.sigma_Nc_off = 0;
        Counts.sigma_Nc_on = 0;
    
        Options.t_step = 1;
        %[HSRLf] = HSRL_retrieval_20220909(Counts,Atmosphere,Options,Spectrum);
        [HSRLf] = backscatterRetrievalMPD(Counts, Model, Spectrum, Options);
        %[HSRLf] = HSRL_retrieval_20230115(Counts,Atmosphere,Options);

        %%%HSRLf.BSRmu = normrnd(HSRLf.BSR,HSRLf.sigma_BR);
        %%%HSRL.fBSR(:,:,jjjj) = fillmissing(HSRLf.BSRmu,'linear');
        HSRL.fBSR(:,:,jjjj)=HSRLf.BSR;

        Counts.Nc_on = Counts.gon(:,:,jjjj);
        Counts.Nc_off = Counts.goff(:,:,jjjj);
        Counts.Nm_on = Counts.gon_mol(:,:,jjjj);
        Counts.Nm_off = Counts.goff_mol(:,:,jjjj);
        Counts.sigma_Nm_off = 0;
        Counts.sigma_Nm_on = 0;
        Counts.sigma_Nc_off = 0;
        Counts.sigma_Nc_on = 0;
    
        Options.t_step = 1;
        %[HSRLg] = HSRL_retrieval_20220909(Counts,Atmosphere,Options,Spectrum);
        [HSRLg] = backscatterRetrievalMPD(Counts, Model, Spectrum, Options);
        %[HSRLg] = HSRL_retrieval_20230115(Counts,Atmosphere,Options);

        % % HSRLg.BSRmu = normrnd(HSRLg.BSR,HSRLg.sigma_BR);
        % % HSRL.gBSR(:,:,jjjj) = fillmissing(HSRLg.BSRmu,'linear');
        HSRL.gBSR(:,:,jjjj) = HSRLg.BSR;

        k = ones(1,1)./(1*1);
        HSRL.fBSR(:,:,jjjj) = nanconv(HSRL.fBSR(:,:,jjjj),k,'edge','nanout');
        HSRL.gBSR(:,:,jjjj) = nanconv(HSRL.gBSR(:,:,jjjj),k,'edge','nanout');
        HSRL.fBSR828(:,:,jjjj) = HSRLf.BSR*828/770;
        HSRL.gBSR828(:,:,jjjj) = HSRLg.BSR*828/770;
        

    elseif strcmp(Options.MPDname,'03') || strcmp(Options.MPDname,'05')
        Counts.Nc_on = Counts.fon(:,:,jjjj);
        Counts.Nc_off = Counts.foff(:,:,jjjj);
        Counts.Nm_on = Counts.fon_mol(:,:,jjjj);
        Counts.Nm_off = Counts.foff_mol(:,:,jjjj);
        [HSRLf] = backscatterRetrievalMPD(Counts, Model, Spectrum, Options);
        HSRL.fBSR(:,:,jjjj) = HSRLf.BSR;
        HSRL.fBSR828(:,:,jjjj) = HSRLf.BSR*828/770;

        Counts.Nc_on = Counts.gon(:,:,jjjj);
        Counts.Nc_off = Counts.goff(:,:,jjjj);
        Counts.Nm_on = Counts.gon_mol(:,:,jjjj);
        Counts.Nm_off = Counts.goff_mol(:,:,jjjj);
        [HSRLg] = backscatterRetrievalMPD(Counts, Model, Spectrum, Options);
        HSRL.gBSR(:,:,jjjj) = HSRLg.BSR;
        HSRL.gBSR828(:,:,jjjj) = HSRLg.BSR*828/770;

    else

        LidarData.Range = Range.rm;
        LidarData.Time = Time.ts;
        LidarData.OfflineCombinedTotalCounts = Counts.foff;
        LidarData.OfflineMolecularTotalCounts = Counts.foff_mol;
        WeatherData.Temperature = Model.T;
        WeatherData.Pressure = Model.P;
        if span_days(1) <  datetime(2021,12,1,1,0,0,'TimeZone','UTC')
            [LidarData]=BackscatterRetrievalBoulder062021(LidarData,WeatherData);
        else
            [LidarData]=BackscatterRetrievalBoulder031122(LidarData,WeatherData);
        end
        HSRL.fBSR(:,:,jjjj) = LidarData.UnmaskedBackscatterRatio;
        HSRL.fBSR(:,:,jjjj) = fillmissing(LidarData.UnmaskedBackscatterRatio,'linear');
        HSRL.fBa = LidarData.UnmaskedAerosolBackscatterCoefficient;
        HSRL.fBm = LidarData.MolecularBackscatterCoefficient;
        %Ba828 = LidarData.UnmaskedAerosolBackscatterCoefficient828;
        %Bm828 = LidarData.MolecularBackscatterCoefficient828;
        HSRL.fBSR828 =LidarData.UnmaskedBackscatterRatio828;
    
        LidarData.OfflineCombinedTotalCounts = Counts.goff;
        LidarData.OfflineMolecularTotalCounts = Counts.goff_mol;
        if span_days(1) <  datetime(2021,12,1,1,0,0,'TimeZone','UTC')
            [LidarData]=BackscatterRetrievalBoulder062021(LidarData,WeatherData);
        else
            [LidarData]=BackscatterRetrievalBoulder031122(LidarData,WeatherData);
        end
        HSRL.gBSR(:,:,jjjj) = LidarData.UnmaskedBackscatterRatio;
        HSRL.gBSR(:,:,jjjj) = fillmissing(LidarData.UnmaskedBackscatterRatio,'linear');
        HSRL.gBa = LidarData.UnmaskedAerosolBackscatterCoefficient;
        HSRL.gBm = LidarData.MolecularBackscatterCoefficient;
        %Ba828 = LidarData.UnmaskedAerosolBackscatterCoefficient828;
        %Bm828 = LidarData.MolecularBackscatterCoefficient828;
        HSRL.gBSR828 =LidarData.UnmaskedBackscatterRatio828;
        HSRL.gBSR(:,:,jjjj) = fillmissing(LidarData.UnmaskedBackscatterRatio828,'linear');
    end
%%

    k = ones(1,1)./(1*1);
    
    %=== apply mask
    
    %==== Smooth HSRL
    HSRL.BSR = nanconv(HSRL.BSR(:,:,end),k,'edge','nanout');
    HSRLfull.BSR = nanconv(HSRLfull.BSR(:,:,end),k,'edge','nanout');
    HSRL.gBSR(:,:,jjjj) = nanconv(HSRL.gBSR(:,:,jjjj),k,'edge','nanout');
    HSRL.fBSR(:,:,jjjj) = nanconv(HSRL.fBSR(:,:,jjjj),k,'edge','nanout');
    
        
%%
    %==========================
    %= Pertabative absorption =
    %==========================
    disp('Calculating absorption')
            
    % === Zeroth Order ===    
    %2nd order error
    [alpha_0_full] = alpha_0(Counts.o2on_noise,Counts.o2off_noise,Range.rangeBin);
    alpha_0_full = real(alpha_0_full);
    
    %Calculate modeled offline absorption
    alpha_0_off = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_offline,Model.WV,Constants);
    Alpha.alpha_0_full = alpha_0_full+alpha_0_off;
    
    %Set mask to model
    nanAlpha = isnan(Alpha.alpha_0_full);
    Alpha.alpha_0_full(nanAlpha) = Model.absorption(nanAlpha);

  
    %Bootsrapping f and g alpha
    Alpha.alpha_0f(:,:,jjjj) = alpha_0(Counts.fon(:,:,jjjj),Counts.foff(:,:,jjjj),Range.rangeBin);
    Alpha.alpha_0f(:,:,jjjj) = real(Alpha.alpha_0f(:,:,jjjj))+alpha_0_off;
    %Alpha.alpha_0f(nanAlpha) = Model.absorption(nanAlpha);
    Alpha.alpha_0g(:,:,jjjj) = alpha_0(Counts.gon(:,:,jjjj),Counts.goff(:,:,jjjj),Range.rangeBin);
    Alpha.alpha_0g(:,:,jjjj) = real(Alpha.alpha_0g(:,:,jjjj))+alpha_0_off;
    %Alpha.alpha_0g(nanAlpha) = Model.absorption(nanAlpha);

    
    % Alpha 0 water vapor
    [Alpha.alpha_0wv] = alpha_0(Counts.wvon_noise,Counts.wvoff_noise,Range.rangeBin);
    Alpha.alpha_0wv = real(Alpha.alpha_0wv);

    %Bootsrapping water vapor alpha
    Alpha.alpha_0wvf(:,:,jjjj) = real(alpha_0(Counts.fwvon(:,:,jjjj),Counts.fwvoff(:,:,jjjj),Range.rangeBin));
    Alpha.alpha_0wvg(:,:,jjjj) = real(alpha_0(Counts.gwvon(:,:,jjjj),Counts.gwvoff(:,:,jjjj),Range.rangeBin));
    
    % Water vapor cross section calculateion
    [~,cross_section,~,~,g_wv] = cross_section_wv_828_model(Model.T,Model.P,Spectrum.nu_wvon,Alpha.alpha_0wv,0,Constants);
    [~,cross_sectionOff,~,~] = cross_section_wv_828_model(Model.T,Model.P,Spectrum.nu_wvoff,Alpha.alpha_0wv,0,Constants);
    % Water vapor number density
    N_wv0 = Alpha.alpha_0wv./(cross_section-cross_sectionOff);

    %Recalculate with water vapor number density
    [~,cross_section,~,~,g_wv] = cross_section_wv_828_model(Model.T,Model.P,Spectrum.nu_wvon,Alpha.alpha_0wv,N_wv0,Constants);
    [~,cross_sectionOff,~,~] = cross_section_wv_828_model(Model.T,Model.P,Spectrum.nu_wvoff,Alpha.alpha_0wv,N_wv0,Constants);
    N_wv0 = Alpha.alpha_0wv./(cross_section-cross_sectionOff);
    
    Alpha.N_wv0=N_wv0;
    
%%
    %==============Set Model WV to DIAL======
    %Purturbative WV
    T_etalonwv_on = HSRL.WVOnlineTransmission./max(HSRL.WVOnlineTransmission);
    %altitude in km
    altitude = 1.5719;
    [Alpha.alpha_1wv, Alpha.alpha_2wv,Spectrum] = pertAbsorptionwv(Alpha.alpha_0wv, T_etalonwv_on, Model, Range, Time, Spectrum, HSRL.BSR828, Options, Constants, altitude);
    
    N_wv = (Alpha.alpha_0wv+Alpha.alpha_1wv+ Alpha.alpha_2wv)./(cross_section-cross_sectionOff);
    [Alpha.alpha_1wvf(:,:,jjjj), Alpha.alpha_2wvf,Spectrum] = pertAbsorptionwv(Alpha.alpha_0wvf(:,:,jjjj), T_etalonwv_on, Model, Range, Time, Spectrum, HSRL.fBSR828(:,:,jjjj), Options, Constants, altitude);
    
    N_wvf(:,:,jjjj) = (Alpha.alpha_0wvf(:,:,jjjj)+Alpha.alpha_1wvf(:,:,jjjj)+ Alpha.alpha_2wvf)./(cross_section-cross_sectionOff);
   
    [Alpha.alpha_1wvg(:,:,jjjj), Alpha.alpha_2wvg,Spectrum] = pertAbsorptionwv(Alpha.alpha_0wvg(:,:,jjjj), T_etalonwv_on, Model, Range, Time, Spectrum, HSRL.gBSR828(:,:,jjjj), Options, Constants, altitude);
    
    N_wvg(:,:,jjjj) = (Alpha.alpha_0wvg(:,:,jjjj)+Alpha.alpha_1wvg(:,:,jjjj)+ Alpha.alpha_2wvg)./(cross_section-cross_sectionOff);
   
    N_wv0f(:,:,jjjj) = Alpha.alpha_0wvf(:,:,jjjj)./(cross_section-cross_sectionOff);
    N_wv0g(:,:,jjjj) = Alpha.alpha_0wvg(:,:,jjjj)./(cross_section-cross_sectionOff);
    
    %Smooth WV
    k = ones(4,6)./(4*6);
    N_wvm = nanconv(N_wv,k,'edge','nanout');
    %N_wvm(cloud_SDm_above)=nan;
    N_wv0m = nanconv(N_wv0,k,'edge','nanout');
    %N_wv0m(cloud_SDm_above)=nan;
    
    AbsHumm = N_wvm.*Constants.mWV*1000; %[g/m3]
    AbsHum0m = N_wv0m.*Constants.mWV*1000; %[g/m3]
    
    AbsHumRawm = N_wv.*Constants.mWV*1000; %[g/m3]
    AbsHumRawm(cloud_SDm_above)=nan;
    AbsHum0Rawm = N_wv0.*Constants.mWV*1000; %[g/m3]
    AbsHum0Rawm(cloud_SDm_above)=nan;
    %%%%%%SET MODEL TO WV RETRIEVAL
    %Model.WV = fillmissing(N_wvm,'linear');
    Model.WV(~isnan(N_wv0m)) = N_wv0m(~isnan(N_wv0m));
    Model.WV = fillmissing(Model.WV,'linear',1);
    Model.WV = fillmissing(Model.WV,'linear',2);
    
    Alpha.AbsHum0Rawm = AbsHum0Rawm;
    Alpha.AbsHumRawm = AbsHumRawm;
    Alpha.N_wv0m=N_wv0m;
    Alpha.N_wvm=N_wvm;
    
    %%
    % === Total alpha WV===    
    Alpha.alpha_total_rawwv = Alpha.alpha_0wv + Alpha.alpha_1wv + Alpha.alpha_2wv;
    N_wv = Alpha.alpha_total_rawwv./cross_section; %[molecule/m3] wv number density
    
    %%
    % === Smoothing zeroth order
    Alpha.alpha_0_full=real(Alpha.alpha_0_full);
    
    T_etalon_on = HSRL.onlineCombinedTransmission./max(HSRL.onlineCombinedTransmission);
    T_etalon_off = HSRL.offlineCombinedTransmission./max(HSRL.offlineCombinedTransmission);
    
%%
    % === Purtabative O2 absorption ===
    [Alpha.alpha_total_rawFull,Alpha.alpha_1,Alpha.alpha_2,~] = pertAbsorption(Alpha.alpha_0_full, T_etalon_on, T_etalon_off, Model, Range, Time, Spectrum, HSRLfull.BSR, Options,true,true,Constants);
  
    Alpha.alpha_total_raw = Alpha.alpha_total_rawFull;
    [Alpha.alpha_total_rawf(:,:,jjjj)] = pertAbsorption(Alpha.alpha_0f(:,:,jjjj), T_etalon_on, T_etalon_off, Model, Range, Time, Spectrum, HSRL.fBSR(:,:,jjjj), Options,true,true,Constants);
    [Alpha.alpha_total_rawg(:,:,jjjj)] = pertAbsorption(Alpha.alpha_0g(:,:,jjjj), T_etalon_on, T_etalon_off, Model, Range, Time, Spectrum, HSRL.gBSR(:,:,jjjj), Options,true,true,Constants);
    
    Alpha.alpha_total_err = zeros(size(Alpha.alpha_total_rawFull));
    
%%
    
    %== Force total alpha to its modeled surface value ==
    %[~,cut] = min(abs(Range.rm-500));             % Index where rm is closest to chosen value
    cut=2;
    Alpha.alpha_total_cut = [Model.absorption(1,:); NaN((cut - 2),Time.i_time); Alpha.alpha_total_rawFull(cut:end,:)];
    Alpha.alpha_total_cut = fillmissing(Alpha.alpha_total_cut,'linear');
    
    %Alpha.alpha_total_cut = Alpha.alpha_0;
    
%%
    %===== Soothing alpha =====
    k=ones(1,1)./(1*1);% Kernel
    
    %Appy cloud mask before smoothing
    Alpha.alpha_total_cut(cloud_SDm_above) = NaN;          % Replace mask with NaNs
    Alpha.alpha_total_cut(cloud_SDm_above) = NaN; 
    
    %Convolve with kernel
    Alpha.alpha_total_filt = nanconv(Alpha.alpha_total_cut,k,'edge','nanout');
    % Alpha.alpha_total_filt(1:UpperSmoothingCut,:) = nanconv(Alpha.alpha_total_cut(1:UpperSmoothingCut,:),k,'edge','nanout');
    % Alpha.alpha_total_filt(UpperSmoothingCut+1:end,:) = nanconv(Alpha.alpha_total_cut(UpperSmoothingCut+1:end,:),upperK,'edge','nanout');

    Alpha.alpha_totals = Alpha.alpha_total_filt;
    
    Alpha.alpha_total_rawFull = nanconv(Alpha.alpha_total_rawFull,k,'edge','nanout');
    Alpha.alpha_total_rawf(:,:,jjjj) = nanconv(Alpha.alpha_total_rawf(:,:,jjjj),k,'edge','nanout');
    Alpha.alpha_total_rawg(:,:,jjjj) = nanconv(Alpha.alpha_total_rawg(:,:,jjjj),k,'edge','nanout');
    Alpha.alpha_total_rawFull = nanconv(Alpha.alpha_total_rawFull,k,'edge','nanout');
    
    Alpha.alpha_0m = Alpha.alpha_0_full;
    Alpha.alpha_0m(cloud_SDm_above)=nan;
    Alpha.alpha_0s= nanconv(Alpha.alpha_0m,k,'edge','nanout');
    
    
    % % k1 = round((size(k,1)-1)./2);
    % % k2 = round((size(k,2)-1)./2);
    % % for iii = 1:size(Alpha.alpha_total_err,1)
    % %     for jjj = 1:size(Alpha.alpha_total_err,2)
    % %         if iii <= k1 || jjj <= k2
    % %             Alpha.alpha_total_errs(iii,jjj) = Alpha.alpha_total_err(iii,jjj);
    % %         elseif iii >= (size(Alpha.alpha_total_err,1)-k1) || jjj >= (size(Alpha.alpha_total_err,2)-k2)
    % %             Alpha.alpha_total_errs(iii,jjj) = Alpha.alpha_total_err(iii,jjj);
    % %         else
    % %             Alpha.alpha_total_errs(iii,jjj) = sqrt(sum(sum((Alpha.alpha_total_err(iii-k1:iii+k1,jjj-k2:jjj+k2)).^2)))./numel(k);
    % %         end
    % %     end
    % % end
    
    
%%
    % apply SNR mask again
    Alpha.alpha_0m = Alpha.alpha_0_full;
    Alpha.alpha_0m(cloud_SDm_above) = NaN;                  % Replace mask with NaNs
    
    Alpha.alpha_total = real(Alpha.alpha_totals);
    Alpha.alpha_totalm = Alpha.alpha_total;
    Alpha.alpha_totalm(cloud_SDm_above) = NaN;          % Replace mask with NaNs
    
%%
    %==========================
    %= Temperature Function =
    %==========================
    disp('Temp retrieval function')
    
    % Testing different models
    %Model.T = Model.Ts + lapseRate .* Range.rm;                           %[K] (1 x r) Temperature model as a function of r 
    %  Model.T = 240*ones(Range.i_range,Time.i_time);
    %  Model.Ts = Model.T(1,:);
    % Model.P = Model.Ps .* (Model.Ts./Model.T).^(-5.2199);                       %[atm] (1 x r) Pressure model as a function of r  
    
    startLapse = Model.lapseRate;
    [Temperature.T_final_testFull(:,:,jjjj),Temperature.L_fit_sm_testFull,Temperature.Ts_fitFull,Temperature.Patm_finalFull,Temperature.mean_lapse_rateFull,Temperature.exclusionFull,Temperature.deltaTFull] =  temperatureRetrieval(Model.T,Time.ts,Range.rm,Model.P,Model.WV,Spectrum.nu_online,Alpha.alpha_total_rawFull,0,cloud_SDm_above|SNRm,Model.Ts,Model.Ps,startLapse,Constants);
    Temperature.L_fit_sm_test = Temperature.L_fit_sm_testFull;
    Temperature.Ts_fit = Temperature.Ts_fitFull;
    Temperature.Patm_final = Temperature.Patm_finalFull;
    Temperature.mean_lapse_rate = Temperature.mean_lapse_rateFull;
    Temperature.exclusion = Temperature.exclusionFull;
    Temperature.deltaT = Temperature.deltaTFull;
    Temperature.T_final_test = Temperature.T_final_testFull;
    [Temperature.T_final_testf(:,:,jjjj),~,~,~,~,~,Temperature.deltaTFullf] =  temperatureRetrieval(Model.T,Time.ts,Range.rm,Model.P,Model.WV,Spectrum.nu_online,Alpha.alpha_total_rawf(:,:,jjjj),0,cloud_SDm_above|SNRm,Model.Ts,Model.Ps,startLapse,Constants);
    [Temperature.T_final_testg(:,:,jjjj),~,~,~,~,~,Temperature.deltaTFullg] =  temperatureRetrieval(Model.T,Time.ts,Range.rm,Model.P,Model.WV,Spectrum.nu_online,Alpha.alpha_total_rawg(:,:,jjjj),0,cloud_SDm_above|SNRm,Model.Ts,Model.Ps,startLapse,Constants);
    
end
 %%   

%bootstrapping temperature and wv error
B = size(Temperature.T_final_testf,3);
tempStd = sqrt((1/(2*(B-1))) *sum((Temperature.T_final_testf-Temperature.T_final_testg).^2,3) );

Temperature.T_final_testFull=Temperature.T_final_testFull(:,:,end);
% mask based on temperature convergence
deltaTMask = false(size(Temperature.T_final_testFull(:,:,end)));
deltaTMask(abs(Temperature.deltaTFull(:,:,end))>=2e-5) = true;
Temperature.T_final_testFull(deltaTMask)=nan;

alpha0wvStd = sqrt((1/(2*(B-1))) *sum((Alpha.alpha_0wvf-Alpha.alpha_0wvg).^2,3) );
N_wv0Std = sqrt((1/(2*(B-1))) *sum((N_wv0f-N_wv0g).^2,3) );
N_wvStd = sqrt((1/(2*(B-1))) *sum((N_wvf-N_wvg).^2,3) );


%%
%======================
%==== Calculate Radiosonde temperature and retrieval temperature difference
%==== (not strickly necessary for temperature retrieval)
%======================
if ~isempty(Sonde.sonde_ind)
    
    TPerfect = 296+Range.rm*-6.5/1000;
    PPerfect = 1 * (296./TPerfect).^(-5.2199);
    WVPerfect = zeros(size(TPerfect));
    %absorptionPerfect = absorption_O2_770_model(TPerfect,PPerfect,Spectrum.nu_online(1),WVPerfect,Constants);
    Temperature.T_final_test2=zeros(Range.i_range,size(Sonde.sonde_ind,2));
    abssorptionSonde = zeros(Range.i_range,size(Sonde.sonde_ind,2));
    Pfinal2 = nan(Range.i_range,size(Sonde.sonde_ind,2));
   for iii = 1:size(Sonde.sonde_ind,2)
        %for iii = 7
        sonde_index = iii;
        p_point = Sonde.sonde_ind(:,sonde_index);
        
        k = ones(1,1)./(1);
        abssorptionSonde(:,iii) = [nanconv(Sonde.absorption_sonde{sonde_index}(1:end-1),k,'edge','nanout') ;0];
        
        abssorptionSonde(:,iii) = Sonde.absorption_sonde{sonde_index};
        
        sondeModelT = Sonde.Tsurf(sonde_index)+Range.rm*-6.5/1000;
        sondeModelP = Sonde.Psurf(sonde_index) .* (Sonde.Tsurf(sonde_index)./sondeModelT).^(-5.2558);
        
        [Tfinal,Lapse,Ts_fit,Pfinal,~,~,~] =  temperatureRetrieval(diag(Model.T(:,p_point)),Time.ts(:,p_point(1)),Range.rm,diag(Model.P(:,p_point)),Sonde.WV_sonde(:,sonde_index),Spectrum.nu_online(:,p_point(1)),abssorptionSonde(:,iii),0,diag(cloud_SDm_above(:,p_point)),Sonde.Tsurf(sonde_index),Sonde.Psurf(sonde_index),startLapse,Constants);

        Tfinal(1:7,:) = NaN;
        Pfinal(1:7,:) = NaN;
        Pfinal2(:,iii) = Pfinal;
        Temperature.T_final_test2(:,iii)=Tfinal;
        
    end
    
    tempComparison = ones(Range.i_range,size(Sonde.sonde_ind,2));
    for jj = 1:size(Sonde.sonde_ind,2)
        tempComparison(:,jj) = Temperature.T_final_test2(:,jj)-Sonde.T_sonde(:,jj);
    end
    meanTemp = mean(tempComparison,2,'omitnan');
    stdTemp = std(tempComparison,0,2,'omitnan');
    
    divider = 12;
    figure(49)
    tempProb = zeros(Range.i_range,121);
    for ii = 1:size(Sonde.sonde_ind,2)
        for i = 1:Range.i_range
            for j = -60:60
                tempDiff = Temperature.T_final_test2(i,ii)-Sonde.T_sonde(i,ii);
                tempProb(i,j+61) = tempProb(i,j+61) + double((tempDiff<(j+1)/divider)&&(tempDiff>=j/divider));
            end
        end
    end
    tempProb(tempProb==0)=nan;
    imAlpha=ones(size(tempProb));
    imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
    imagesc((-60:60)/divider,Range.rkm,tempProb,'AlphaData',imAlpha)
    hold on
    plot(meanTemp,Range.rkm,'b','linewidth',2)
    plot(meanTemp+stdTemp,Range.rkm,'--k')
    plot(meanTemp-stdTemp,Range.rkm,'--k')
    hold off
    colormap(flipud(hot))%colormap hot
    set(gca,'Color','#D3D3D3')
    a = colorbar;
    a.Label.String = 'Occurrences';
    hold on
    xline(0)
    hold off
    ylim([0 4])
    xlim([-3 3])
    set(gca, 'YDir','normal')
    xlabel('\DeltaT MPD-Sonde (^oC)')
    ylabel('Range (km)')
    grid on
    title('temperature')
    
    figure(50)
    int = 1;
    hist11 = zeros(size(Sonde.sonde_ind,2)*size(Sonde.sonde_ind,1));
    for ii = 1:size(Sonde.sonde_ind,2)
        for jj = 1:size(Sonde.sonde_ind,1)
            hist11(int)=Temperature.T_final_test2(jj,ii)-Sonde.T_sonde(jj,ii);
            int=int+1;
        end
    end
    meanHist = mean(hist11,'omitnan');
    stdHist = std(hist11,'omitnan');
    histogram(hist11,'BinWidth',.01)
    title(sprintf('Histogram T_{DIAL}-T_{sonde}\n Mean %0.5f, std %0.5f',meanHist(1),stdHist(1)))

    figure(53)
    tiledlayout(1,2)
    nexttile
    plot(Temperature.T_final_test2-Sonde.T_sonde,Range.rkm)
    hold on
    %plot(meanTemp,Range.rm,'linewidth',3)
    %title(sprintf('temperature\n mean %0.3d std %0.3d',meanHist(1),stdHist(1)))
    %title(sprintf('Temperature Difference'))
    xlabel('T_{retrieval}-T_{sonde} (K)')
    xlabel('\Delta T (retrieval-sonde) (^oC)','FontWeight','bold')
    ylabel('Range (km)','FontWeight','bold')
    %xlim([-10 10])
    %xlim([-.01  .04])
    hold off
    grid on

   % figure(54)
    nexttile
    plot(Pfinal2-Sonde.P_sonde,Range.rkm)
    %title('Pressure Difference')
    xlabel('\Delta P (retrieval-sonde) (atm)','FontWeight','bold')
    %ylabel('Range (km)')
    grid on
    
    figure(5544)
    subplot(2,1,1)
    for ii = 1:size(Sonde.sonde_ind,2)
        plot(Time.thr(Sonde.sonde_ind(1,ii)),Model.Ts(Sonde.sonde_ind(1,ii))- Sonde.Tsurf(ii),'.')
        hold on
    end
    ylabel('surface weather - sonde (K)')
    hold off
    subplot(2,1,2)
    for ii = 1:size(Sonde.sonde_ind,2)
        plot(Time.thr(Sonde.sonde_ind(1,ii)),Model.Ps(Sonde.sonde_ind(1,ii))- Sonde.Psurf(ii),'.')
        hold on
    end
    ylabel('surface weather - sonde (atm)')
    xlabel('Time (hr)')
    hold off
end
%%

%=== Cut temperature to surface value ====
% Temperature.T_final_test_cut = [Model.Ts(1,:); NaN((cut - 2),Time.i_time); Temperature.T_final_test(cut:end,:)];
% Temperature.T_final_test_cut = fillmissing(Temperature.T_final_test_cut,'linear');

%Appy cloud mask before smoothing
%Temperature.T_final_test_cut(cloud_SDm_above) = NaN;          % Replace mask with NaNs

%Temperature.T_final_testf(isnan(gg)) = NaN;
%Temperature.T_final_testg(isnan(gg)) = NaN;
%[Ez,Et,minSigz,minSigt] = findMinE(Temperature.T_final_testf,Temperature.T_final_testg,0);
%[Temperature.T_final_tests2] = applyFilter(minSigz,minSigt,Temperature.T_final_test_cut);


%====================
%=== apply mask
%==== Smooth temperature
%====================
k = ones(4,6)./(4*6);
Temperature.T_final_testFulls = nanconv(Temperature.T_final_testFull(:,:,end),k,'edge','nanout');
Temperature.T_final_testfs = zeros(size(Temperature.T_final_testFull,1),size(Temperature.T_final_testFull,2),iter);
Temperature.T_final_testgs = zeros(size(Temperature.T_final_testFull,1),size(Temperature.T_final_testFull,2),iter);
for iii = 1:size(Temperature.T_final_testf,3)   
    Temperature.T_final_testfs(:,:,iii) = nanconv(Temperature.T_final_testf(:,:,iii),k,'edge','nanout');
    Temperature.T_final_testgs(:,:,iii) = nanconv(Temperature.T_final_testg(:,:,iii),k,'edge','nanout');
end

Temperature.TempStd = std(Temperature.T_final_testf(:,:,:),0,3);
k1 = round((size(k,1)-1)./2);
k2 = round((size(k,2)-1)./2);
tempStdss = nan(size(Temperature.TempStd));
tempStdssS = nan(size(Temperature.TempStd));
tempStdssSl = nan(size(Temperature.TempStd));
for iii = 1:size(Temperature.TempStd,1)
    for jjj = 1:size(Temperature.TempStd,2)
        if iii <= k1 || jjj <=k2
        Temperature.TempStds(iii,jjj) = Temperature.TempStd(iii,jjj);
            tempStdssS(iii,jjj) = tempStd(iii,jjj,end);
            tempStdssSl(iii,jjj) = tempStd(iii,jjj,end);
        elseif iii >= (size(Temperature.TempStd,1)-k1) || jjj >= (size(Temperature.TempStd,2)-k2)
            Temperature.TempStds(iii,jjj) = Temperature.TempStd(iii,jjj);

            tempStdss(iii,jjj) = tempStd(iii,jjj,end);
            tempStdssS(iii,jjj) = tempStd(iii,jjj,end);
            tempStdssSl(iii,jjj) = tempStd(iii,jjj,end);
        else
            Temperature.TempStds(iii,jjj) = sqrt(sumsqr(Temperature.TempStd(iii-k1:iii+k1,jjj-k2:jjj+k2)))./numel(k);

            tempStdss(iii,jjj) = sqrt(sumsqr(tempStd(iii-k1:iii+k1,jjj-k2:jjj+k2,end)))./numel(k);
            tempStdssS(iii,jjj) = sqrt(sumsqr(tempStd(iii-k1:iii+k1-1,jjj-k2:jjj+k2-1,end)))./numel(k);
            tempStdssSl(iii,jjj) = sqrt(sumsqr(tempStd(iii-k1+1:iii+k1,jjj-k2+1:jjj+k2,end)))./numel(k);

            tempStdssSl(iii,jjj) = sqrt(sumsqr(tempStd(iii-k1+1:iii+k1,jjj-k2+1:jjj+k2,end)))./numel(k);  
        end    
    end
end

tempStdssSS = sqrt(nanconv(tempStd(:,:,end).^2,k.*numel(k),'edge','nanout').*numel(k))./numel(k);

%tempStds = sqrt((1/(2*(B-1))) *sum((Temperature.T_final_testfs-Temperature.T_final_testgs).^2,3) );

tempStd = sqrt((1./(2.*(permute(1:B,[1 3 2])-1))) .*cumsum((Temperature.T_final_testf-Temperature.T_final_testg).^2,3) );

tempStds = tempStdssSS;
Temperature.TempStds = tempStds;

%%
%=====Smooth watervapor and calculate absolute humidity
k = ones(2,1)./(2*1);
N_wv0s = nanconv(N_wv0(:,:,end),k.*numel(k),'edge','nanout');
N_wvs = nanconv(N_wv(:,:,end),k.*numel(k),'edge','nanout');
N_wvStds = sqrt(nanconv(N_wvStd(:,:,end).^2,k.*numel(k),'edge','nanout').*numel(k))./numel(k);
N_wv0Stds = sqrt(nanconv(N_wv0Std(:,:,end).^2,k.*numel(k),'edge','nanout').*numel(k))./numel(k);

AbsHum0s = N_wv0s.*Constants.mWV*1000; %[g/m^3]
AbsHums = N_wvs.*Constants.mWV*1000;
AbsHum0Stds = N_wv0Stds.*Constants.mWV*1000;
AbsHumStds = N_wvStds.*Constants.mWV*1000;
%%
%=== apply mask
%Temperature.T_finalm = Temperature.T_final_tests ;
%Temperature.T_finalm(cloud_SDm_above) = NaN;

Temperature.T_finalm = Temperature.T_final_testFulls;

%%

Results.Temperature = Temperature.T_finalm;
Results.Range = Range.rm;
Results.Time = Time.ts;
Results.Date = Time.date_ts;
Results.BSR = HSRL.BSR;
Results.WV = N_wv;
Results.O2absorption_0 = Alpha.alpha_0m;
Results.O2absorption_1 = Alpha.alpha_1;
Results.O2absorption_2 = Alpha.alpha_2;
Results.O2absorption_total = Alpha.alpha_totalm;
Results.TemperatureError = Temperature.TempStds;
% Results.TemperatureErrorMask = mask;
% Results.AbsoluteHumidity = AbsHum0m;
% Results.UnmaskedBSR = HSRLfull.BSR;

%%
%===============
%=== Figures ===
%===============

%=Tick spacing and formating
tickHours =12;
tickMin = 0;
[y,m,d]=ymd(span_days(1));
date2=datetime(y,m,d,tickHours,tickMin,0);
date1=datetime(y,m,d,0,0,0);
Format.tickSpacing = datenum(date2)-datenum(date1);
Format.tickAngle = 30;
Format.dateTickFormat ='mm/dd HH:MM';
%Format.dateTickFormat ='mm/dd';

%= Plot time for profiles
%plot_time = datetime(2023,2,13,23,0,0,'TimeZone','UTC');%yyyy,mm,dd,hh,mm
%plot_time = datetime(2023,8,16,18,0,0,'TimeZone','UTC');%yyyy,mm,dd,hh,mm
%plot_time = datetime(2023,8,5,17,0,0,'TimeZone','UTC');%yyyy,mm,dd,hh,mm
plot_time = datetime(2024,7,19,12,00,0,'TimeZone','UTC');%yyyy,mm,dd,hh,mm
%plot_time = datetime(2022,6,2,16,00,0,'TimeZone','UTC');%yyyy,mm,dd,hh,mm
[~,p_point] = min(abs(plot_time-Time.date_ts)); % Find closest value to 338min for comparison to other program
p_point(1:length(Range.rm),1)=p_point;

%= Plot time for profiles with sondes
sonde_index = 1;
%p_point = Sonde.sonde_ind(:,sonde_index);

%mask = logical(Temperature.TempStds>2) | cloud_SDm_above;
%mask = logical(tempStds>2) | cloud_SDm_above;
mask = logical(tempStds>5) | isnan(tempStds);
wvmask = logical(AbsHumStds>5) | isnan(AbsHumStds);

%mask = cloud_SDm_above;
Temperature.T_finalm(mask) = nan;
Alpha.Alpha_totalm(mask)=nan;
AbsHum0sm = AbsHum0s;
AbsHum0sm(wvmask)=nan;

[~,lowAltindex] = min(abs(Range.rm-lowAlt));
AbsHum0sm(1:lowAltindex,:)=nan;
AbsHumsm = AbsHums;
AbsHumsm(wvmask)=nan;
AbsHumsm(1:lowAltindex,:)=nan;
%HSRL.BSR(mask)=nan;

Alpha.alpha_totals(mask) = nan;

plot_O2(p_point,sonde_index,span_days,Sonde,Model,Counts,Range,Time,Options,Temperature,Format,Alpha,cloud_SDm,HSRL,Data,cloud_SDm,mask,N_wv,N_wv0,N_wvm,N_wv0m,AbsHumm,AbsHum0m,Constants)