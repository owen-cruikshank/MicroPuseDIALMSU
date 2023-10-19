% O2_DIAL_Temp_retrieval_odcV2.m
% Analysis program for O2 DIAL instrument from 
clear all

%=Date Range

date_start = datetime(2022,4,21,'TimeZone','UTC');%yyyy,mm,dd
date_end = datetime(2022,7,16,'TimeZone','UTC');%yyyy,mm,dd

date_start = datetime(2022,6,22,'TimeZone','UTC');%yyyy,mm,dd
date_end = datetime(2022,6,22,'TimeZone','UTC');%yyyy,mm,dd

date_start = datetime(2023,8,1,'TimeZone','UTC');%yyyy,mm,dd
date_end = datetime(2023,8,1,'TimeZone','UTC');%yyyy,mm,dd

span_days = date_start:date_end;

%=Time and range averaging
Options.intTime = 1;  %[min] Integration time
Options.intTime = 10;  %[min] Integration time
Options.intRange = 1; %[bins] Integration range

Options.t_avg = 1;     %[bins] Time smoothing bins
Options.oversample = 1; %[bins] Range smoothing bins

%====================
%==== Constants =====
%====================
Constant.g0 = 9.80665;       %[m/s^2] Gravitational acceleration 
Constant.M_air = 0.0289644;  %[kg/mol] Molar mass of Earth's air 
Constant.R = 8.3144598;      %[J/(mol*K)] Universal gas constant 

Constant.c = 2.99792458E8;           %[m/s] Speed of light 
Constant.kb = 1.38065E-23;           %[J/K][m^2 kg s-2 K-1] Boltzman's constant 
Constant.h = 6.626E-34;              %[Js] Planck's constant 
Constant.mo2 = 5.314E-26;            %[kg] Mass O2 molecule 
Constant.mWV = 2.9915e-26;           %[kg] Mass H2O molecule
Constant.m_air = 4.792E-26;          %[kg] Mass of air
Constant.q_O2 = .2095;               %[unitless] O2 atmospheric mixing ratio
Constant.No = 2.47937E25;            %[1/m^3] Loschmidt's number  (referenced to 296 K and 1 atm)

%======================
%==== Load data
%==== photon counts, radiosondes, weather station data, HSRL data, etc.
%======================

%==Load data from MSU instument==
%[Range,Time,Counts,Sonde,Model,Spectrum,HSRL,Data] = loadMSUdata(span_days,Options,Constant);

% ==Load data from SGP ARM instument==
%date_begin = datetime(2019,4,17); date_end   = datetime(2019,4,22);
%span_days = date_begin:date_end;        % Days in set [datetime]
%[Range,Time,Counts,Sonde,Model,Spectrum,BSR,Data] = loadSGPdata(span_days,Options,Constant);

%==Load data from Boulder instument==

%Options.MPDname = '05';
%Options.MPDname = '01';
Options.MPDname = '03';
Options.DataPath = 'C:\Users\Owen\OneDrive - Montana State University\Research\O2 DIAL\Data';
[Range,Time,Counts,Sonde,Model,Spectrum,HSRL,Data] =loadBoulderdata3(span_days,Options,Constant);


%%
%======================
%= Cloud and SNR mask =
%======================
disp('Calculating masks')
 %(hours) Time to plot mask data
cloud_p_point = 20;

%Threshold for initial SNR mask
SNR_threshold = 2300;
%SNR_threshold = 1000;

%Threshold for Cloud mask
%SD_threshold = 0.5;  
SD_threshold = 1.5;
SD_threshold = 5;

BGmult =1; % Multiplier for background for SNR calculation
%Low altitude masking
lowAlt = 330;
[SNRm , cloud_SDm_above, cloud_SDm,~] = mask_O2_BSR(cloud_p_point,SNR_threshold,SD_threshold,Options.oversample,Options.t_avg,Counts,BGmult,Time,Range,HSRL.BSR,lowAlt);

cloud_SDm_above(cloud_SDm_above==-1)=0;

cloud_SDm_above = ~(logical(cloud_SDm_above));
SNRm = ~(logical(SNRm));

cloud_SDm_above = cloud_SDm_above | SNRm;
cloud_SDm = logical(cloud_SDm);
%clear  o2on_SNR SNRm
%%
Counts.o2on(cloud_SDm_above) = nan;
Counts.o2off(cloud_SDm_above) = nan;
Counts.o2on_mol(cloud_SDm_above) = nan;
Counts.o2off_mol(cloud_SDm_above) = nan;
Counts.wvon(cloud_SDm_above) = nan;
Counts.wvoff(cloud_SDm_above) = nan;
%%
%Number of poisson thinning iterations
iter = 20;
Counts = poissonThin2(Counts,cloud_SDm_above,iter);

%%
for jjjj = 1:iter
    if jjjj >=2
%         Model.T = real(fillmissing(Temperature.T_finalm,'linear'));
%         Model.P = real(fillmissing(Temperature.Patm_final,'linear'));

        Model.T = fillmissing(Temperature.L_fit_sm_test(:,:,end).*Range.rm+Temperature.Ts_fit(:,:,end),'linear');
        Model.Ts =Temperature.Ts_fit(:,:,end);

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
        Pg = fillmissing(Model.Ps.*(Temperature.Ts_fit(:,:,end)./(Temperature.Ts_fit(:,:,end)+Temperature.L_fit_sm_test(:,:,end).*Range.rm)).^(gamma./Temperature.L_fit_sm_test(:,:,end)),'linear');
      
        %Pg = fillmissing(Model.Ps.*(Model.Ts./(Model.Ts+LapseRand.*Range.rkm)).^(gamma./LapseRand/1000),'linear');
        Model.P = Pg;

        Model.absorption = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_online,Model.WV);

        Model.absorption_off = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_offline,Model.WV); %[m-1] Funcrtion to calculate theoretical absorption
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
        [HSRL] = HSRL_retrieval_20220909(Counts,Atmosphere,Options);
        %[HSRL] = HSRL_retrieval_20230115(Counts,Atmosphere,Options);
    
    elseif strcmp(Options.MPDname,'03')
        Counts.Nc_on = Counts.o2on;%.*Counts.NBins;
        Counts.Nc_off = Counts.o2off;%.*Counts.NBins;
        Counts.Nm_on = Counts.o2on_mol;%.*Counts.NBins;
        Counts.Nm_off = Counts.o2off_mol;%.*Counts.NBins;
        [HSRL] = backscatterRetrievalMPD03(Counts, Model, Spectrum);
        % Counts.Nc_on = Counts.o2on_bgsub;%.*Counts.NBins;
        % Counts.Nc_off = Counts.o2off_bgsub;%.*Counts.NBins;
        % Counts.Nm_on = Counts.o2on_bgsub_mol;%.*Counts.NBins;
        % Counts.Nm_off = Counts.o2off_bgsub_mol;%.*Counts.NBins;
        % [HSRLfull] = backscatterRetrievalMPD03(Counts, Model, Spectrum);
    
    end
    
    %Calulate HSRL at 828nm
    HSRL.Bm828 = HSRL.Bm *(770/828)^4;
    HSRL.Ba828 = HSRL.Ba*(770/828);
    HSRL.BSR828 = HSRL.Ba828./HSRL.Bm828+1;
  
%%
    if strcmp(Options.MPDname,'MSU')

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
        %[HSRLf] = HSRL_retrieval_20220909(Counts,Atmosphere,Options);
        [HSRLf] = HSRL_retrieval_20230115(Counts,Atmosphere,Options);

        HSRLf.BSRmu = normrnd(HSRLf.BSR,HSRLf.sigma_BR);
        HSRL.fBSR(:,:,jjjj) = fillmissing(HSRLf.BSRmu,'linear');

        Counts.Nc_on = Counts.gon(:,:,jjjj);
        Counts.Nc_off = Counts.goff(:,:,jjjj);
        Counts.Nm_on = Counts.gon_mol(:,:,jjjj);
        Counts.Nm_off = Counts.goff_mol(:,:,jjjj);
        Counts.sigma_Nm_off = 0;
        Counts.sigma_Nm_on = 0;
        Counts.sigma_Nc_off = 0;
        Counts.sigma_Nc_on = 0;
    
        Options.t_step = 1;
        %[HSRLg] = HSRL_retrieval_20220909(Counts,Atmosphere,Options);
        [HSRLg] = HSRL_retrieval_20230115(Counts,Atmosphere,Options);

        HSRLg.BSRmu = normrnd(HSRLg.BSR,HSRLg.sigma_BR);
        HSRL.gBSR(:,:,jjjj) = fillmissing(HSRLg.BSRmu,'linear');

    elseif strcmp(Options.MPDname,'03')
        Counts.Nc_on = Counts.fon(:,:,jjjj);
        Counts.Nc_off = Counts.foff(:,:,jjjj);
        Counts.Nm_on = Counts.fon_mol(:,:,jjjj);
        Counts.Nm_off = Counts.foff_mol(:,:,jjjj);
        [HSRLf] = backscatterRetrievalMPD03(Counts, Model, Spectrum);
        HSRL.fBSR(:,:,jjjj) = HSRLf.BSR;

        Counts.Nc_on = Counts.gon(:,:,jjjj);
        Counts.Nc_off = Counts.goff(:,:,jjjj);
        Counts.Nm_on = Counts.gon_mol(:,:,jjjj);
        Counts.Nm_off = Counts.goff_mol(:,:,jjjj);
        [HSRLg] = backscatterRetrievalMPD03(Counts, Model, Spectrum);
        HSRL.gBSR(:,:,jjjj) = HSRLg.BSR;

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
    % ==== Afterpulse correction ====
    % pulseON = 4.146e+07*10*(Range.rm).^(-2.536*1);
    % % pulseOFF = 4.146e+07*1*(rm_over).^(-2.536*1.0);
    % pulseOFF = 3.682e+07*10*1.3*(Range.rm).^(-2.512*1);
    % 
    % % pulseON = 4.146e+07*10*1.3*(Range.rm).^(-2.536*1);
    % % % pulseOFF = 4.146e+07*1*(rm_over).^(-2.536*1.0);
    % % pulseOFF = 3.682e+07*10*1*(Range.rm).^(-2.512*1);
    % Counts.o2on = Counts.o2on-pulseON;
    % Counts.o2off = Counts.o2off-pulseOFF;
    
%%
    %==========================
    %= Pertabative absorption =
    %==========================
    disp('Calculating absorption')
            
    % === Zeroth Order ===
    ind_r_lo = 1:Range.i_range-Options.oversample;                                            % High range vector
    ind_r_hi = 1+Options.oversample:Range.i_range;                                            % Low range vector
    ln_o2 = log((Counts.o2on(ind_r_lo,:).*Counts.o2off(ind_r_hi,:))./(Counts.o2on(ind_r_hi,:).*Counts.o2off(ind_r_lo,:))); % Natural log of counts
        
    Alpha.alpha_0_raw = ln_o2./2./(Range.rangeBin*Options.oversample);                              %[1/m] 
    %Alpha.alpha_0 = interp2(Time.ts,Range.rm(ind_r_lo),Alpha.alpha_0_raw,Time.ts,Range.rm);
    Alpha.alpha_0 = interp2(Time.ts,Range.rm(ind_r_lo)+Range.rangeBin./2,Alpha.alpha_0_raw,Time.ts,Range.rm);
    Alpha.alpha_0_raw = real(Alpha.alpha_0_raw);
    Alpha.alpha_0 = real(Alpha.alpha_0);
    alpha_0 = fillmissing(Alpha.alpha_0,'nearest');
    
    %2nd order error
        alpha_0 = zeros(size(Counts.o2on));
        alpha_0(1,:,:) = log((Counts.o2on(1,:).*Counts.o2off(2,:))./(Counts.o2on(2,:).*Counts.o2off(1,:)))./2./Range.rangeBin;
        for iii = 2:length(Range.rm)-1
            %dg1(iii,:,:) = (g1(iii+1,:,:)-g1(iii-1,:,:))/2/Range.rangeBin;
            alpha_0(iii,:)=log((Counts.o2on(iii-1,:).*Counts.o2off(iii+1,:))./(Counts.o2on(iii+1,:).*Counts.o2off(iii-1,:)))./2./Range.rangeBin./2;
        end
        alpha_0(end,:,:) = log((Counts.o2on(end-1,:).*Counts.o2off(end,:))./(Counts.o2on(end,:).*Counts.o2off(end-1,:)))./2./Range.rangeBin;
    
        alpha_0 = real(alpha_0);
    
    alpha_0_off = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_offline,Model.WV);
    
    Alpha.alpha_0 = alpha_0-alpha_0_off;
    
    nanAlpha = isnan(Alpha.alpha_0);
    Alpha.alpha_0(nanAlpha) = Model.absorption(nanAlpha);
    
    Counts.o2off_mol_cor = Counts.o2off_mol.*HSRL.BSR;
        alpha_0_mol = zeros(size(Counts.o2on));
        alpha_0_mol(1,:,:) = log((Counts.o2on_mol(1,:).*Counts.o2off_mol_cor(2,:))./(Counts.o2on_mol(2,:).*Counts.o2off_mol_cor(1,:)))./2./Range.rangeBin;
        for iii = 2:length(Range.rm)-1
            alpha_0_mol(iii,:)=log((Counts.o2on_mol(iii-1,:).*Counts.o2off_mol_cor(iii+1,:))./(Counts.o2on_mol(iii+1,:).*Counts.o2off_mol_cor(iii-1,:)))./2./Range.rangeBin./2;
        end
        alpha_0_mol(end,:,:) = log((Counts.o2on_mol(end-1,:).*Counts.o2off_mol_cor(end,:))./(Counts.o2on_mol(end,:).*Counts.o2off_mol_cor(end-1,:)))./2./Range.rangeBin;
    
    ln_o2 = log((Counts.fon(ind_r_lo,:,jjjj).*Counts.foff(ind_r_hi,:,jjjj))./(Counts.fon(ind_r_hi,:,jjjj).*Counts.foff(ind_r_lo,:,jjjj))); % Natural log of counts
        
    Alpha.alpha_0f_raw = real(ln_o2)./2./(Range.rangeBin*Options.oversample);                              %[1/m] 
    %Alpha.alpha_0 = interp2(Time.ts,Range.rm(ind_r_lo),Alpha.alpha_0_raw,Time.ts,Range.rm);
    Alpha.alpha_0f(:,:,jjjj) = interp2(Time.ts,Range.rm(ind_r_lo)+Range.rangeBin./2,Alpha.alpha_0f_raw,Time.ts,Range.rm);
    Alpha.alpha_0f(:,:,jjjj) = fillmissing(Alpha.alpha_0f(:,:,jjjj),'nearest');
    
    ln_o2 = log((Counts.gon(ind_r_lo,:,jjjj).*Counts.goff(ind_r_hi,:,jjjj))./(Counts.gon(ind_r_hi,:,jjjj).*Counts.goff(ind_r_lo,:,jjjj))); % Natural log of counts
        
    Alpha.alpha_0g_raw = real(ln_o2)./2./(Range.rangeBin*Options.oversample);                              %[1/m] 
    %Alpha.alpha_0 = interp2(Time.ts,Range.rm(ind_r_lo),Alpha.alpha_0_raw,Time.ts,Range.rm);
    Alpha.alpha_0g(:,:,jjjj) = interp2(Time.ts,Range.rm(ind_r_lo)+Range.rangeBin./2,Alpha.alpha_0g_raw,Time.ts,Range.rm);
    Alpha.alpha_0g(:,:,jjjj) = fillmissing(Alpha.alpha_0g(:,:,jjjj),'nearest');
    
    
    %central diff
    % don = zeros(size(Counts.o2on));
    % don(1,:) = log(Counts.o2on(2,:))-log(Counts.o2on(1,:));
    % doff = zeros(size(Counts.o2on));
    % doff(1,:) = log(Counts.o2off(2,:))-log(Counts.o2off(1,:));
    % for iii = 2:length(Range.rm)-1
    % don(iii,:) = (log(Counts.o2on(iii+1,:))-log(Counts.o2on(iii-1,:)))./2;
    % doff(iii,:) = (log(Counts.o2off(iii+1,:))-log(Counts.o2off(iii-1,:)))./2;
    % end
    % don(end,:) = log(Counts.o2on(end,:))-log(Counts.o2on(end-1,:));
    % doff(end,:) = log(Counts.o2off(end,:))-log(Counts.o2off(end-1,:));
    % 
    % alpha_0_center = (doff-don)./2./Range.rm;
    %Alpha.alpha_0 = alpha_0_center;
    
    %%%%%Alpha.alpha_0 = Alpha.alpha_0 - Model.absorption_off; %Subtract offline absorption
    
    
    %zero order error alpha
    Alpha.alpha_0_err = (1./2./(Range.rangeBin*Options.oversample)) .* sqrt(1./Counts.NBins) .*  sqrt(...
     (-sqrt(Counts.o2on(ind_r_lo,:)+Counts.bg_o2on)./Counts.o2on(ind_r_lo,:)).^2 ...
    +(sqrt(Counts.o2on(ind_r_hi,:)+Counts.bg_o2on)./Counts.o2on(ind_r_hi,:)).^2 ...
    +(sqrt(Counts.o2off(ind_r_lo,:)+Counts.bg_o2off)./Counts.o2off(ind_r_lo,:)).^2 ...
    +(-sqrt(Counts.o2off(ind_r_hi,:)+Counts.bg_o2off)./Counts.o2off(ind_r_hi,:)).^2);
    
    Alpha.alpha_0_err(end+1,:,:)=Alpha.alpha_0_err(end,:,:);
    
    
    ln_o2 = log((Counts.wvon(ind_r_lo,:).*Counts.wvoff(ind_r_hi,:))./(Counts.wvon(ind_r_hi,:).*Counts.wvoff(ind_r_lo,:))); % Natural log of counts  
    Alpha.alpha_0wv_raw = ln_o2./2./(Range.rangeBin*Options.oversample); 
    %Alpha.alpha_0wv = interp2(Time.ts,Range.rm(ind_r_lo),Alpha.alpha_0wv_raw,Time.ts,Range.rm);
    Alpha.alpha_0wv = interp2(Time.ts,Range.rm(ind_r_lo)+Range.rangeBin./2,Alpha.alpha_0wv_raw,Time.ts,Range.rm);
    Alpha.alpha_0wv = fillmissing(Alpha.alpha_0wv,'nearest');
    
    clear ln_o2 foff fon goff gon
    
    [~,cross_section,~,~] = cross_section_wv_828_model(Model.T,Model.P,Spectrum.nu_wvon,Alpha.alpha_0wv);
    
    [~,cross_sectionOff,~,~] = cross_section_wv_828_model(Model.T,Model.P,Spectrum.nu_wvoff,Alpha.alpha_0wv);
    
    N_wv0 = Alpha.alpha_0wv./(cross_section-cross_sectionOff);
    
    Alpha.N_wv0=N_wv0;
    
%%
    %==============Set Model WV to DIAL======
    
    %Purturbative WV
    % load(fullfile('CalibrationData','TransmittanceData.mat'))
    % T_etalon_on = double(interp1(double(OnlineWavelength)*10^9,OnlineCombinedTransmittance,Spectrum.lambda_scan_3D_short));
    % T_etalon_off = double(interp1(double(OfflineWavelength)*10^9,OfflineCombinedTransmittance,Spectrum.lambda_scan_3D_short_off));
    % 
    T_etalon_on = HSRL.onlineCombinedTransmission./max(HSRL.onlineCombinedTransmission);


    %altitude in km
    altitude = 1.5719;
    
    [Alpha.alpha_1wv, Alpha.alpha_2wv,Spectrum] = pertAbsorptionwv(Alpha.alpha_0wv, T_etalon_on, Model, Range, Time, Spectrum, HSRL.BSR828, ind_r_lo,ind_r_hi, Options, Constant, altitude);
    
    N_wv = (Alpha.alpha_0wv+Alpha.alpha_1wv+ Alpha.alpha_2wv)./(cross_section-cross_sectionOff);
    
    %Smooth WV
    % k = ones(4,8)./(4*8);     % Kernel
    % %k = ones(3,7)./(3*7);     % Kernel
    % k = ones(2,2)./(2*2);     % Kernel
    k = ones(3,3)./(3*3);     % Kernel
    
    
    N_wvm = nanconv(N_wv,k,'edge','nanout');
    N_wvm(cloud_SDm_above)=nan;
    
    N_wv0m = nanconv(N_wv0,k,'edge','nanout');
    N_wv0m(cloud_SDm_above)=nan;
    
    AbsHumm = N_wvm.*Constant.mWV*1000; %[g/m3]
    AbsHum0m = N_wv0m.*Constant.mWV*1000; %[g/m3]
    
    AbsHumRawm = N_wv.*Constant.mWV*1000; %[g/m3]
    AbsHumRawm(cloud_SDm_above)=nan;
    AbsHum0Rawm = N_wv0.*Constant.mWV*1000; %[g/m3]
    AbsHum0Rawm(cloud_SDm_above)=nan;
    %%%%%%SET MODEL TO WV RETRIEVAL
    Model.WV = fillmissing(N_wvm,'linear');
    Model.WV = fillmissing(N_wv0m,'linear');
    
    Alpha.AbsHum0Rawm = AbsHum0Rawm;
    Alpha.AbsHumRawm = AbsHumRawm;
    Alpha.N_wv0m=N_wv0m;
    Alpha.N_wvm=N_wvm;
    
    %%
    % === Total alpha ===
    
    Alpha.alpha_total_rawwv = Alpha.alpha_0wv + Alpha.alpha_1wv + Alpha.alpha_2wv;
    N_wv = Alpha.alpha_total_rawwv./cross_section; %[molecule/m3] wv number density
    
    % === Smoothing zeroth order
    
    Alpha.alpha_0=real(Alpha.alpha_0);
    
    % === Molecular alpha calcultion =====
    
%%
    % === Purtabative absorption ===
    [Alpha.alpha_total_raw,Alpha.alpha_1,Alpha.alpha_2,Spectrum] = pertAbsorption(Alpha.alpha_0, T_etalon_on, Model, Range, Time, Spectrum, HSRL.BSR, ind_r_lo,ind_r_hi, Options,true);
    
    [Alpha.alpha_total_rawf(:,:,jjjj)] = pertAbsorption(Alpha.alpha_0f(:,:,jjjj), T_etalon_on, Model, Range, Time, Spectrum, HSRL.fBSR(:,:,jjjj), ind_r_lo,ind_r_hi, Options,true);
    [Alpha.alpha_total_rawg(:,:,jjjj)] = pertAbsorption(Alpha.alpha_0g(:,:,jjjj), T_etalon_on, Model, Range, Time, Spectrum, HSRL.gBSR(:,:,jjjj), ind_r_lo,ind_r_hi, Options,true);
    
%%
    Alpha.alpha_total_err = zeros(size(Alpha.alpha_total_raw));
    
%%
    
    %== Force total alpha to its modeled surface value ==
    %[~,cut] = min(abs(Range.rm-500));             % Index where rm is closest to chosen value
    cut=2;
    Alpha.alpha_total_cut = [Model.absorption(1,:); NaN((cut - 2),Time.i_time); Alpha.alpha_total_raw(cut:end,:)];
    Alpha.alpha_total_cut = fillmissing(Alpha.alpha_total_cut,'linear');
    
    %Alpha.alpha_total_cut = Alpha.alpha_0;
    
%%
    %===== Soothing alpha =====
    k = ones(2,2)./(2*2);% Kernel
    
    %Appy cloud mask before smoothing
    Alpha.alpha_total_cut(cloud_SDm_above) = NaN;          % Replace mask with NaNs
    Alpha.alpha_total_cut(cloud_SDm_above) = NaN; 
    % Alpha.alpha_total_rawmf = Alpha.alpha_total_rawf;
    % Alpha.alpha_total_rawmg = Alpha.alpha_total_rawg;
    % Alpha.alpha_total_rawmf(isnan(gg))=nan;
    % Alpha.alpha_total_rawmg(isnan(gg))=nan;
    
    [~,UpperSmoothingCut] = min(abs(Range.rm-3000));
    % upperK  = ones(7,21)./(7*21);
    % upperK = ones(1,1)./(1);
    upperK = k;
    
    %Convolve with kernel
    Alpha.alpha_total_filt = nanconv(Alpha.alpha_total_cut,k,'edge','nanout');
    Alpha.alpha_total_filt(1:UpperSmoothingCut,:) = nanconv(Alpha.alpha_total_cut(1:UpperSmoothingCut,:),k,'edge','nanout');
    Alpha.alpha_total_filt(UpperSmoothingCut+1:end,:) = nanconv(Alpha.alpha_total_cut(UpperSmoothingCut+1:end,:),upperK,'edge','nanout');
    Alpha.alpha_totals = Alpha.alpha_total_filt;
    
    Alpha.alpha_total_rawf(:,:,jjjj) = nanconv(Alpha.alpha_total_rawf(:,:,jjjj),k,'edge','nanout');
    Alpha.alpha_total_rawg(:,:,jjjj) = nanconv(Alpha.alpha_total_rawg(:,:,jjjj),k,'edge','nanout');
    
    Alpha.alpha_0m = Alpha.alpha_0;
    Alpha.alpha_0m(cloud_SDm_above)=nan;
    Alpha.alpha_0s= nanconv(Alpha.alpha_0m,k,'edge','nanout');
    
    
    k1 = round((size(k,1)-1)./2);
    k2 = round((size(k,2)-1)./2);
    for iii = 1:size(Alpha.alpha_total_err,1)
        for jjj = 1:size(Alpha.alpha_total_err,2)
            if iii <= k1 || jjj <= k2
                Alpha.alpha_total_errs(iii,jjj) = Alpha.alpha_total_err(iii,jjj);
            elseif iii >= (size(Alpha.alpha_total_err,1)-k1) || jjj >= (size(Alpha.alpha_total_err,2)-k2)
                Alpha.alpha_total_errs(iii,jjj) = Alpha.alpha_total_err(iii,jjj);
            else
                %Alpha.alpha_total_errs(iii,jjj) = sqrt(sumsqr(Alpha.alpha_total_err(iii-k1:iii+k1,jjj-k2:jjj+k2)))./numel(k);
                Alpha.alpha_total_errs(iii,jjj) = sqrt(sum(sum((Alpha.alpha_total_err(iii-k1:iii+k1,jjj-k2:jjj+k2)).^2)))./numel(k);
                %Alpha.alpha_total_errs(iii,jjj) = sqrt(sum(sum((Alpha.alpha_total_err(iii-k1:iii+k1,jjj-k2:jjj+k2)).^2))./numel(k));
            end
        end
    end
    
    % Alpha.alpha_total_filt = nanconv(Alpha.alpha_0,k,'edge','nanout');
    % Alpha.alpha_totals = Alpha.alpha_total_filt;
    
    %%%Alpha.alpha_totals = Alpha.alpha_total_filt2;
    
    % Alpha.alpha_total_filtf = nanconv(Alpha.alpha_total_rawf,k,'edge','nanout');
    % Alpha.alpha_totalsf = Alpha.alpha_total_filtf;
    % Alpha.alpha_total_filtg = nanconv(Alpha.alpha_total_rawg,k,'edge','nanout');
    % Alpha.alpha_totalsg = Alpha.alpha_total_filtg;
    
    
%%
    % apply SNR mask again
    Alpha.alpha_0m = Alpha.alpha_0;
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
    [Temperature.T_final_test(:,:,jjjj),Temperature.L_fit_sm_test,Temperature.Ts_fit,Temperature.Patm_final,Temperature.mean_lapse_rate,Temperature.exclusion,Temperature.Titer] =  temperatureRetrieval(Model.T,Time.ts,Range.rm,Model.P,Model.WV,Spectrum.nu_online,Alpha.alpha_totals,0,cloud_SDm_above|SNRm,Model.Ts,Model.Ps,startLapse);
    
    [Temperature.T_final_testf(:,:,jjjj)] =  temperatureRetrieval(Model.T,Time.ts,Range.rm,Model.P,Model.WV,Spectrum.nu_online,Alpha.alpha_total_rawf(:,:,jjjj),0,cloud_SDm_above|SNRm,Model.Ts,Model.Ps,startLapse);
    [Temperature.T_final_testg(:,:,jjjj)] =  temperatureRetrieval(Model.T,Time.ts,Range.rm,Model.P,Model.WV,Spectrum.nu_online,Alpha.alpha_total_rawg(:,:,jjjj),0,cloud_SDm_above|SNRm,Model.Ts,Model.Ps,startLapse);
end

%Temperature.T_final_test = (mean(Temperature.T_final_testf,3)+mean(Temperature.T_final_testg,3))/2;

%bootstrapping 
B = size(Temperature.T_final_testf,3);
tempStd = sqrt((1/(2*(B-1))) *sum((Temperature.T_final_testf-Temperature.T_final_testg).^2,3) );

%%%%Temperature.TempSTD = 0;
%%

Temperature.T_final_test0 =  temperatureRetrieval(Model.T,Time.ts,Range.rm,Model.P,Model.WV,Spectrum.nu_online,Alpha.alpha_0s,0,cloud_SDm_above,Model.Ts,Model.Ps,startLapse);
%%
if ~isempty(Sonde.sonde_ind)
    
    TPerfect = 296+Range.rm*-6.5/1000;
    PPerfect = 1 * (296./TPerfect).^(-5.2199);
    WVPerfect = zeros(size(TPerfect));
    absorptionPerfect = absorption_O2_770_model(TPerfect,PPerfect,Spectrum.nu_online(1),WVPerfect);
    Temperature.T_final_test2=zeros(Range.i_range,size(Sonde.sonde_ind,2));
    abssorptionSonde = zeros(Range.i_range,size(Sonde.sonde_ind,2));
    Pfinal2 = zeros(Range.i_range,size(Sonde.sonde_ind,2));
    for iii = 1:size(Sonde.sonde_ind,2)
        sonde_index = iii;
        p_point = Sonde.sonde_ind(:,sonde_index);
        
        %k = ones(4,1)./(4);
        k = ones(1,1)./(1);
        abssorptionSonde(:,iii) = [nanconv(Sonde.absorption_sonde{sonde_index}(1:end-1),k,'edge','nanout') ;0];
        
        %abssorptionSonde(:,iii) = filter2(k,Sonde.absorption_sonde{sonde_index},'same');
        %abssorptionSonde = Sonde.absorption_sonde{sonde_index};
        %abssorptionSonde = absorptionPerfect;
        
        sondeModelT = Sonde.Tsurf(sonde_index)+Range.rm*-6.5/1000;
        sondeModelP = Sonde.Psurf(sonde_index) .* (Sonde.Tsurf(sonde_index)./sondeModelT).^(-5.2558);
        
        %[Tfinal,~,~,Pfinal,~,~,~] =  temperatureRetrieval(Model.T(:,p_point(1)),Time.ts(:,p_point(1)),Range.rm,Model.P(:,p_point(1)),Sonde.WV_sonde(:,sonde_index),Spectrum.nu_online(:,p_point(1)),abssorptionSonde,SNRm(:,p_point(1)),cloud_SDm_above(:,p_point(1)),Model.Ts(:,p_point(1)),Model.Ps(:,p_point(1)));
        %%[Tfinal,~,~,Pfinal,~,~,~] =  temperatureRetrieval(diag(Model.T(:,p_point)),Time.ts(:,p_point(1)),Range.rm,diag(Model.P(:,p_point)),Sonde.WV_sonde(:,sonde_index),Spectrum.nu_online(:,p_point(1)),abssorptionSonde,diag(SNRm(:,p_point)),diag(cloud_SDm_above(:,p_point)),Sonde.Tsurf(sonde_index),Sonde.Psurf(sonde_index));
        [Tfinal,Lapse,Ts_fit,Pfinal,~,~,~] =  temperatureRetrieval(diag(Model.T(:,p_point)),Time.ts(:,p_point(1)),Range.rm,diag(Model.P(:,p_point)),Sonde.WV_sonde(:,sonde_index),Spectrum.nu_online(:,p_point(1)),abssorptionSonde(:,iii),0,diag(cloud_SDm_above(:,p_point)),Sonde.Tsurf(sonde_index),Sonde.Psurf(sonde_index),startLapse);
        %[Tfinal,Lapse,Ts_fit,Pfinal,~,~,Titer(:,iii,:)] =  temperatureRetrieval(diag(Model.T(:,p_point)),Time.ts(:,p_point(1)),Range.rm,diag(Model.P(:,p_point)),Sonde.WV_sonde(:,sonde_index),Spectrum.nu_online(:,p_point(1)),abssorptionSonde(:,iii),0,diag(cloud_SDm_above(:,p_point)),Model.Ts(Sonde.sonde_ind(1,sonde_index)),(Model.Ps(Sonde.sonde_ind(1,sonde_index))));
        %%%[Tfinal,Lapse,Ts_fit,Pfinal,~,~,Titer(:,iii,:)] =  temperatureRetrieval(sondeModelT,Time.ts(:,p_point(1)),Range.rm,sondeModelP,Sonde.WV_sonde(:,sonde_index),Spectrum.nu_online(:,p_point(1)),abssorptionSonde,diag(SNRm(:,p_point)),diag(cloud_SDm_above(:,p_point)),Sonde.Tsurf(sonde_index),Sonde.Psurf(sonde_index));
        %%[Tfinal,Lapse,Ts_fit,Pfinal,~,~,~] =  temperatureRetrieval(sondeModelT,Time.ts(:,p_point(1)),Range.rm,sondeModelP,zeros(size(Sonde.WV_sonde(:,sonde_index))),Spectrum.nu_online(:,p_point(1)),abssorptionSonde,diag(SNRm(:,p_point)),diag(cloud_SDm_above(:,p_point)),Sonde.T_sonde(1,sonde_index),Sonde.P_sonde(1,sonde_index));
        %%[Tfinal,Lapse,Ts_fit,Pfinal,~,~,~] =  temperatureRetrieval(sondeModelT,Time.ts(:,p_point(1)),Range.rm,sondeModelP,diag(N_wvm(:,Sonde.sonde_ind(:,sonde_index))),Spectrum.nu_online(:,p_point(1)),abssorptionSonde,diag(SNRm(:,p_point)),diag(cloud_SDm_above(:,p_point)),Sonde.T_sonde(1,sonde_index),Sonde.P_sonde(1,sonde_index));
        %[Tfinal,~,~,Pfinal,~,~,~] =  temperatureRetrieval(sondeModelT,Time.ts(:,p_point(1)),Range.rm,sondeModelP,WVPerfect,Spectrum.nu_online(:,p_point(1)),abssorptionSonde,diag(SNRm(:,p_point)),diag(cloud_SDm_above(:,p_point)),Sonde.Tsurf(sonde_index),Sonde.Psurf(sonde_index));
        %[Tfinal,~,~,Pfinal,~,~,~] =  temperatureRetrieval(Model.T(:,p_point(1)),Time.ts(:,p_point(1)),Range.rm,Model.P(:,p_point(1)),WVPerfect,Spectrum.nu_online(1),absorptionPerfect,SNRm(:,p_point(1)),cloud_SDm_above(:,p_point(1)),296,1+.001);
        
        Tfinal(diag(cloud_SDm_above(:,p_point))) = NaN;
        Pfinal(diag(cloud_SDm_above(:,p_point))) = NaN;
        Pfinal2(:,iii) = Pfinal;
        Temperature.T_final_test2(:,iii)=Tfinal;
        
        
        %%%[Temperature.T_final_testf,Temperature.L_fit_sm_test,Temperature.Ts_fit,Temperature.Patm_final,Temperature.mean_lapse_rate,Temperature.exclusion,Temperature.Titer] =  temperatureRetrieval(Model.T,Time.ts,Range.rm,Model.P,Model.WV,Spectrum.nu_online,Alpha.alpha_totalsf,SNRm,cloud_SDm_above);
        %%%[Temperature.T_final_testg,Temperature.L_fit_sm_test,Temperature.Ts_fit,Temperature.Patm_final,Temperature.mean_lapse_rate,Temperature.exclusion,Temperature.Titer] =  temperatureRetrieval(Model.T,Time.ts,Range.rm,Model.P,Model.WV,Spectrum.nu_online,Alpha.alpha_totalsg,SNRm,cloud_SDm_above);
    
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
    plot(Temperature.T_final_test2-Sonde.T_sonde,Range.rkm)
    hold on
    %plot(meanTemp,Range.rm,'linewidth',3)
    title(sprintf('temperature\n mean %0.3d std %0.3d',meanHist(1),stdHist(1)))
    title(sprintf('Temperature Difference'))
    xlabel('T_{retrieval}-T_{sonde} (K)')
    ylabel('Range (km)')
    xlim([-10 10])
    hold off
    grid on

    figure(54)
    plot(Pfinal2-Sonde.P_sonde,Range.rkm)
    title('Pressure Difference')
    xlabel('P_{retrieval}-P_{sonde} (atm)')
    ylabel('Range (km)')
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
    
    
    
    % sonde_index=1;
    % figure(51)
    % plot(Temperature.T_final_test2(:,sonde_index),Range.rkm)
    % hold on
    % plot(permute(Titer(:,sonde_index,:),[1 3 2]),Range.rkm)
    % plot(Sonde.T_sonde(:,sonde_index),Range.rkm,'--')
    % plot(Sonde.Tsurf(1,sonde_index),0,'*')
    % plot(Model.Ts(1,Sonde.sonde_ind(1,sonde_index)),0,'+')
    % title('Temparature')
    % hold off
    % ylim([0 5])
    % 
    % figure(511)
    % plot(Pfinal2(:,sonde_index),Range.rkm)
    % hold on
    % plot(Sonde.P_sonde(:,sonde_index),Range.rkm,'--')
    % plot(Sonde.Psurf(1,sonde_index),0,'*')
    % plot(Model.Ps(1,Sonde.sonde_ind(1,sonde_index)),0,'+')
    % plot(Model.P(:,Sonde.sonde_ind(1,sonde_index)),Range.rkm,'+')
    % hold off
    % title('Pressure')
    % ylim([0 5])
    
    % figure(5111)
    % plot(Sonde.absorption_sonde{sonde_index},Range.rkm)
    % hold on
    % plot(abssorptionSonde(:,sonde_index),Range.rkm)
    % hold off
    % legend('unsmoothed','smoothed')
    % title('absorption')
    % ylim([0 5])
    
    % figure(52)
    % plot(Temperature.T_final_test2(:,sonde_index)-Sonde.T_sonde(:,sonde_index),Range.rkm)
    % xlim([-10 10])
    % ylim([0 5])
    % title('temperature')
    % xline(1)
    % xline(-1)
    % 
    % figure(522)
    % plot(Pfinal2(:,sonde_index)-Sonde.P_sonde(:,sonde_index),Range.rkm)
    % title('Pressure')
    % %xlim([-10 10])
    % ylim([0 5])
    % %xline(1)
    % %xline(-1)
    
    
    
    % figure
    % plot(Pfinal,Range.rkm)
    % hold on
    % plot(Sonde.P_sonde(:,sonde_index),Range.rkm)
    % hold off
    % 
    % figure
    % plot(TPerfect,Range.rm)
    % hold on
    % plot(Tfinal,Range.rm)
    % hold off
    % 
    % figure
    % plot(TPerfect-Tfinal,Range.rm)
    % hold on
    % 
    % hold off

end
%%

%=== Cut temperature to surface value ====
% Temperature.T_final_test_cut = [Model.Ts(1,:); NaN((cut - 2),Time.i_time); Temperature.T_final_test(cut:end,:)];
% Temperature.T_final_test_cut = fillmissing(Temperature.T_final_test_cut,'linear');

%Appy cloud mask before smoothing
Temperature.T_final_test_cut(cloud_SDm_above) = NaN;          % Replace mask with NaNs


%Temperature.T_final_testf(isnan(gg)) = NaN;
%Temperature.T_final_testg(isnan(gg)) = NaN;
%[Ez,Et,minSigz,minSigt] = findMinE(Temperature.T_final_testf,Temperature.T_final_testg,0);
%[Temperature.T_final_tests2] = applyFilter(minSigz,minSigt,Temperature.T_final_test_cut);


k = ones(3,3)./(3*3);     % Kernel

k = ones(4,6)./(4*6);

%=== apply mask

%==== Smooth temperature
Temperature.T_final_tests = nanconv(Temperature.T_final_test(:,:,end),k,'edge','nanout');

Temperature.T_final_testfs = zeros(size(Temperature.T_final_tests,1),size(Temperature.T_final_tests,2),iter);
Temperature.T_final_testgs = zeros(size(Temperature.T_final_tests,1),size(Temperature.T_final_tests,2),iter);
for iii = 1:size(Temperature.T_final_testf,3)
Temperature.T_final_testfs(:,:,iii) = nanconv(Temperature.T_final_testf(:,:,1),k,'edge','nanout');
Temperature.T_final_testgs(:,:,iii) = nanconv(Temperature.T_final_testg(:,:,1),k,'edge','nanout');
end

Temperature.TempStd = zeros(size(Temperature.T_final_tests));
Temperature.TempStd = std(Temperature.T_final_testf(:,:,:),0,3);
%Temperature.TempStds = Temperature.TempStd./21;
k1 = round((size(k,1)-1)./2);
k2 = round((size(k,2)-1)./2);
tempStdss = nan(size(Temperature.TempStd));
for iii = 1:size(Temperature.TempStd,1)
    for jjj = 1:size(Temperature.TempStd,2)

        if iii <= k1 || jjj <=k2
        Temperature.TempStds(iii,jjj) = Temperature.TempStd(iii,jjj);
        elseif iii >= (size(Temperature.TempStd,1)-k1) || jjj >= (size(Temperature.TempStd,2)-k2)
            Temperature.TempStds(iii,jjj) = Temperature.TempStd(iii,jjj);

            tempStdss(iii,jjj) = tempStd(iii,jjj);
        else
            Temperature.TempStds(iii,jjj) = sqrt(sumsqr(Temperature.TempStd(iii-k1:iii+k1,jjj-k2:jjj+k2)))./numel(k);

            tempStdss(iii,jjj) = sqrt(sumsqr(tempStd(iii-k1:iii+k1,jjj-k2:jjj+k2)))./numel(k);
        end
        
    end
end


tempStds = sqrt((1/(2*(B-1))) *sum((Temperature.T_final_testfs-Temperature.T_final_testgs).^2,3) );


tempStd = sqrt((1./(2.*(permute(1:B,[1 3 2])-1))) .*cumsum((Temperature.T_final_testf-Temperature.T_final_testg).^2,3) );

tempStdDiff = diff(tempStd,1,3);

%%
%=== apply mask
Temperature.T_finalm(:,:,jjjj) = Temperature.T_final_tests ;
%Temperature.T_finalm(cloud_SDm_above,end) = NaN;
Temperature.T_finalm2 =Temperature.T_finalm;
Temperature.T_finalm = Temperature.T_finalm(:,:,jjjj);
Temperature.T_finalm(cloud_SDm_above) = NaN;

Temperature.T_final_tests0 = Temperature.T_final_test0;
Temperature.T_final_tests0(cloud_SDm_above) = nan;
Temperature.T_final_tests0 = nanconv(Temperature.T_final_tests0,k,'edge','nanout');

% %=== apply mask
% Temperature.T_finalm = Temperature.T_final_tests ;
% Temperature.T_finalm(cloud_SDm_above) = NaN;



%%

% % % % for m = 1:9
% % % %     for n = 1:4
% % % %          k = ones(n,m)./(n*m);     % Kernel
% % % % 
% % % % 
% % % % 
% % % % %==== Smooth temperature
% % % % Temperature.T_final_tests = nanconv(Temperature.T_final_test(:,:,1),k,'edge','nanout');
% % % % 
% % % % 
% % % % 
% % % % Temperature.TempStd = zeros(size(Temperature.T_final_tests));
% % % % Temperature.TempStd = std(Temperature.T_final_testf(:,:,:),0,3);
% % % % %Temperature.TempStds = Temperature.TempStd./21;
% % % % k1 = round((size(k,1)-1)./2);
% % % % k2 = round((size(k,2)-1)./2);
% % % % for iii = 1:size(Temperature.TempStd,1)
% % % %     for jjj = 1:size(Temperature.TempStd,2)
% % % % 
% % % %         if iii <= k1 || jjj <=k2
% % % %         Temperature.TempStds(iii,jjj) = Temperature.TempStd(iii,jjj);
% % % %         elseif iii >= (size(Temperature.TempStd,1)-k1) || jjj >= (size(Temperature.TempStd,2)-k2)
% % % %             Temperature.TempStds(iii,jjj) = Temperature.TempStd(iii,jjj);
% % % %         else
% % % %             Temperature.TempStds(iii,jjj) = sqrt(sumsqr(Temperature.TempStd(iii-k1:iii+k1,jjj-k2:jjj+k2)))./numel(k);
% % % %         end
% % % %         
% % % %     end
% % % % end
% % % % 
% % % % 
% % % % Temperature.T_finalm = Temperature.T_final_tests ;
% % % % Temperature.T_finalm(cloud_SDm_above) = NaN;
% % % % 
% % % % tempComparison = nan(Range.i_range,size(Sonde.sonde_ind,2));
% % % % for jj = 1:size(Sonde.sonde_ind,2)
% % % %     if jj==12
% % % %     else
% % % %     tempComparison(:,jj) = (diag(Temperature.T_finalm(:,Sonde.sonde_ind(:,jj)))-Sonde.T_sonde(:,jj))./Sonde.T_sonde(:,jj)*100;
% % % %     end
% % % % end
% % % % if ~isempty(Sonde.sonde_ind)
% % % %  meanTemp = mean(tempComparison,2,'omitnan');
% % % %  stdTemp = std(tempComparison,0,2,'omitnan');
% % % % else
% % % %     meanTemp = nan;
% % % %     stdTemp = nan;
% % % % end
% % % % 
% % % % 
% % % % figure(2345)
% % % % plot(Temperature.TempStd(:,p_point(1)),Range.rkm)
% % % % hold on
% % % % plot(Temperature.TempStds(:,p_point(1)),Range.rkm)
% % % % plot(stdTemp,Range.rkm)
% % % % hold off
% % % % xlabel('errest')
% % % %     ylabel('Range (km)')
% % % %     legend('Temperature std', 'Temp std smooth','MDP-Sonde std')
% % % %     grid on
% % % %     title(sprintf('Temp %dx%d',n,m))
% % % % 
% % % % 
% % % % 
% % % %     saveas(gcf,['C:\Users\Owen\OneDrive - Montana State University\Research\Reports\8_25_22\' sprintf('temp%dx%d.png',n,m)])
% % % %     end
% % % % end

%%
% for iiii = 1:size(Sonde.sonde_ind,2)
% 
% sonde_index = iiii;
% p_point = Sonde.sonde_ind(:,sonde_index);
% 
% figure
% errorbar(diag(Temperature.T_finalm(:,p_point))-273,Range.rkm,diag(Temperature.TempStds(:,p_point)))
% errorbar(Range.rkm,diag(Temperature.T_finalm(:,p_point))-273,diag(Temperature.TempStds(:,p_point)))
% hold on
% plot(Range.rkm,Sonde.T_sonde(:,sonde_index)-273)
% xlabel('Range (km)')
% ylabel('Temperature (^oC)')
% legend('T_{MPD}','T_{radiosonde}')
% title(sprintf(['Temperature\n' datestr(Time.date_ts(p_point(1)))]))
% xlim([0 4])
% view([90 -90])
% grid on
% 
% saveas(gcf,['D:\OneDrive - Montana State University\Research\Reports\9_1_22\' sprintf('temp%d.png',iiii)])
% end


figure(54524)
for iii = 1:size(Sonde.sonde_ind,2)

    sonde_index = iii;
p_point = Sonde.sonde_ind(:,sonde_index);
    plot(diag(Temperature.Patm_final(:,p_point))-Sonde.P_sonde(:,sonde_index),Range.rkm)
    hold on
end

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
plot_time = datetime(2023,2,13,23,0,0,'TimeZone','UTC');%yyyy,mm,dd,hh,mm
plot_time = datetime(2023,8,16,18,0,0,'TimeZone','UTC');%yyyy,mm,dd,hh,mm
plot_time = datetime(2023,8,2,12,0,0,'TimeZone','UTC');%yyyy,mm,dd,hh,mm
[~,p_point] = min(abs(plot_time-Time.date_ts)); % Find closest value to 338min for comparison to other program
p_point(1:length(Range.rm),1)=p_point;

%= Plot time for profiles with sondes
sonde_index = 1;
%p_point = Sonde.sonde_ind(:,sonde_index);

mask = logical(Temperature.TempStds>2) | cloud_SDm_above;

%mask = cloud_SDm_above;
Temperature.T_finalm(mask) = nan;
Alpha.Alpha_totalm(mask)=nan;
%HSRL.BSR(mask)=nan;

Alpha.alpha_totals(mask) = nan;

plot_O2(p_point,sonde_index,span_days,Sonde,Model,Counts,Range,Time,Options,Temperature,Format,Alpha,cloud_SDm,HSRL,Data,cloud_SDm,mask,N_wv,N_wv0,N_wvm,N_wv0m,AbsHumm,AbsHum0m)