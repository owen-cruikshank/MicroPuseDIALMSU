function [Range,Time,Counts,Sonde,Model,Spectrum,HSRL,Data] = loadBoulderdata3(span_days,Options,Constant)

disp('Reading in files')
%Read data from MSU DIAL

% cd ../
% Options.path = fullfile(pwd,'Data','NCAR Boulder Data',['mpd_' Options.MPDname '_data']);
% %pathPython = fullfile(pwd,'Data','NCAR Boulder Data','Python');
% pathSonde = fullfile(pwd,'Data','NCAR Boulder Data','Soundings','Soundings-20211220T155456Z-001','Soundings');
% pathSonde = fullfile(pwd,'Data','NCAR Boulder Data');
% cd( fullfile('analysis'))

Options.path = fullfile(Options.DataPath,'NCAR Boulder Data',['mpd_' Options.MPDname '_data']);
%pathPython = fullfile(pwd,'Data','NCAR Boulder Data','Python');
pathSonde = fullfile(Options.DataPath,'NCAR Boulder Data','Soundings','Soundings-20211220T155456Z-001','Soundings');
pathSonde = fullfile(Options.DataPath,'NCAR Boulder Data');



Options.BinTotal = 490;
[Data, Options] = loadMSUNETcdf(span_days,Options);

%%%load('PrePrecipdata041522.mat','Options','Data','Counts')
%%

% for jj = 1:length(span_days)
%     [DataPython] = loadNCARBoulderDataPython(span_days(jj),pathPython);
%     DataSonde = loadNCARBoulderDataSonde(span_days(jj),pathSonde);
% 
%     T_sgp_raw = {[]};
%     rm_sgp = {[]};
%     P_sgp_raw = {[]};
%     datetime_sgp_cell = {[]};
% 
%     if jj==1
% 
%     range_Backscatter_Ratio = DataPython.range;
%     time_Backscatter_Ratio = DataPython.time';
%     Backscatter_Ratio = DataPython.Backscatter_Ratio;
%     Aerosol_Backscatter_Coefficient = DataPython.Aerosol_Backscatter_Coefficient;
%     Backscatter_Ratio_mask = DataPython.Backscatter_Ratio_mask;
%     HSRLMolecular_RayleighBrillioun = DataPython.HSRLMolecular_RayleighBrillioun;
%     r_freqency = DataPython.r_freqency;
% 
%     Absolute_Humidity = DataPython.Absolute_Humidity;
%     Absolute_Humidity_mask = DataPython.Absolute_Humidity_mask;
%     time_Absolute_Humidity = DataPython.time';
%     range_Absolute_Humidity = DataPython.range;
% 
%     time_Temperature = DataPython.time';
%     range_Temperature = DataPython.range;
%     Temperature = DataPython.Temperature_Model;
%     Surface_Temperature_HSRL = DataPython.Surface_Temperature';
% 
%     time_Pressure = DataPython.time';
%     range_Pressure = DataPython.range;
%     Pressure = DataPython.Pressure_Model;
%     Surface_Pressure_HSRL = DataPython.Surface_Pressure';
% 
%     else
%         range_Backscatter_Ratio = DataPython.range;
%         time_Backscatter_Ratio = [time_Backscatter_Ratio DataPython.time'+time_Backscatter_Ratio(end)];
%         Backscatter_Ratio = [Backscatter_Ratio DataPython.Backscatter_Ratio];
%         Aerosol_Backscatter_Coefficient = [Aerosol_Backscatter_Coefficient DataPython.Aerosol_Backscatter_Coefficient];
%         Backscatter_Ratio_mask = [Backscatter_Ratio_mask DataPython.Backscatter_Ratio_mask];
%         %HSRLMolecular_RayleighBrillioun = [HSRLMolecular_RayleighBrillioun DataPython.HSRLMolecular_RayleighBrillioun];
%         HSRLMolecular_RayleighBrillioun = cat(3,HSRLMolecular_RayleighBrillioun,DataPython.HSRLMolecular_RayleighBrillioun);
%         r_freqency = [r_freqency DataPython.r_freqency];
% 
%         Absolute_Humidity = [Absolute_Humidity DataPython.Absolute_Humidity];
%         Absolute_Humidity_mask = [Absolute_Humidity_mask DataPython.Absolute_Humidity_mask];
%         time_Absolute_Humidity = [time_Absolute_Humidity DataPython.time'+time_Absolute_Humidity(end)];
%         range_Absolute_Humidity = DataPython.range;
% 
%         time_Temperature = [time_Temperature DataPython.time'+time_Temperature(end)];
%         range_Temperature = DataPython.range;
%         Temperature = [Temperature DataPython.Temperature_Model];
%         Surface_Temperature_HSRL = [Surface_Temperature_HSRL DataPython.Surface_Temperature'];
% 
%         time_Pressure = [time_Pressure DataPython.time'+time_Pressure(end)];
%         range_Pressure = DataPython.range;
%         Pressure = [Pressure DataPython.Pressure_Model];
%         Surface_Pressure_HSRL = [Surface_Pressure_HSRL DataPython.Surface_Pressure'];
% 
%         
%     end
% end
% 
% Counts.NBinsPython = ones(size(time_Temperature));
% 
% Spectrum.HSRLMolecular_RayleighBrillioun = HSRLMolecular_RayleighBrillioun;
% Spectrum.r_freqency = r_freqency;
% HSRL.Backscatter_Ratio_mask = Backscatter_Ratio_mask;


%%disp('Creating time and range vectors')
% ====================
% === Time Vectors ===
% ====================
%Time.ts = Options.TimeGrid*60*60;
Time.ts = Data.MCS.Channel0.NewTimeGrid*60*60; 
Time.t_step = Time.ts(2)-Time.ts(1);
Time.i_time = length(Time.ts);
Time.date_ts = span_days(1) + seconds(Time.ts);
Time.date_ts.TimeZone = 'UTC';
Time.thr = Time.ts/60/60;

% =====================
% === Range Vectors ===
% =====================
% Create range vector from length of data bins
Range.nsPerBin = 250; %[ns] bin length in nanosections
%Range.NBins = floor(560/Options.intRange); %number of range bins
Range.NBins = floor(490/Options.intRange); %number of range bins
Range.rangeBin = (Constant.c * Range.nsPerBin(1)*10^-9)/2; %range bin length

Range.rm_raw_o2 = 0:Range.rangeBin:Range.NBins(1)*Range.rangeBin+0-Range.rangeBin;    %[m] Create range vector
Range.rm_raw_o2 = -150:Range.rangeBin:Range.NBins(1)*Range.rangeBin-150-Range.rangeBin;    %[m] Create range vector
Range.rm_raw_o2 = -75:Range.rangeBin:Range.NBins(1)*Range.rangeBin-75-Range.rangeBin;    %[m] Create range vector
Range.rm_raw_o2 = -75-30/2:Range.rangeBin:Range.NBins(1)*Range.rangeBin-75-30/2-Range.rangeBin;    %[m] Create range vector
%Range.rm_raw_o2 = -75-30/2-500:Range.rangeBin:Range.NBins(1)*Range.rangeBin-75-30/2-Range.rangeBin-500;
%Range.rm_raw_o2 = -300:Range.rangeBin:Range.NBins(1)*Range.rangeBin-300-Range.rangeBin;
%%%Range.rm_raw_o2 = -75-75:Range.rangeBin:Range.NBins(1)*Range.rangeBin-75-75-Range.rangeBin;    %[m] Create range vector
Range.rm_raw_o2 = Range.rm_raw_o2(:);                           %[m] Convert range vector to column vector
Range.r_max = 6000;                                       %[m] Max range 
Range.rm = Range.rm_raw_o2(Range.rm_raw_o2<=Range.r_max & Range.rm_raw_o2>0);     %[m] Shorten range vector to max range


Range.rm = Range.rm(1:2:end);%integrate to new
Range.rangeBin = Range.rangeBin*2;
Range.rkm = Range.rm/1000;

Range.i_range = length(Range.rm);                               %[none] Size of range vector
%%

disp('Calculating model')

% Calculating temperature and pressure model
            
% Model.Ts = interp1(time_Temperature,Surface_Temperature_HSRL,Time.ts,'nearest',nan); %[K] Surface temperature             
% Model.Ps = interp1(time_Pressure,Surface_Pressure_HSRL,Time.ts,'nearest',nan); %[atm] Absolute surface pressure
Model.Ts = Data.WS.all.Temperature+273.15;
Model.Ps = Data.WS.all.Pressure*0.000987;

lapseRate = -6.5;                                   %[K/km] Guess adiabatic lapse rate  typically -6.5 up to 10km
%%%lapseRate = -8; 
lapseRate = lapseRate / 1000;                       %[K/m] 
Model.lapseRate = lapseRate;
Model.T = Model.Ts + lapseRate .* Range.rm;                           %[K] (1 x r) Temperature model as a function of r 

% Model.T = interp2(double(time_Temperature),double(range_Temperature),double(Temperature),Time.ts,Range.rm,'nearest'); % Interpolate Temperature model to range and time to match other data
% Model.T = fillmissing(Model.T,'nearest',1);
% Model.T = fillmissing(Model.T,'nearest',2);
%%%%T(1,:) = Ts;


Model.P = Model.Ps .* (Model.Ts./Model.T).^(-5.2199);                       %[atm] (1 x r) Pressure model as a function of r  
g0 = 9.80665;               %[m/s/s]Gravitational acceleration
M_air = 0.0289644;          %[kg/mol] Molar mass of air
kB = 1.38065E-23;           %[J/K] Boltzman's constant 
N_A = 6.02214076e23;        %[1/mol] Avagadro's number
R = kB * N_A;               %[J/K/mol] universal gas constant
gamma = g0 * M_air / R;     %[K/m]gravity molar mass of air and gas constant
Model.P = Model.Ps .* (Model.Ts./Model.T).^(gamma./lapseRate);

% Model.P = interp2(double(time_Pressure),double(range_Pressure),double(Pressure),Time.ts,Range.rm,'nearest'); % Interpolate pressure model to range and time to match other data
% Model.P = fillmissing(Model.P,'nearest',1);
% Model.P = fillmissing(Model.P,'nearest',2);
%%%P(1,:) = Ps;

%===============
%= Water Vapor =
%===============

% Profile
% Absolute_Humidity = Absolute_Humidity.*(1-double(Absolute_Humidity_mask));
% Absolute_Humidity(Absolute_Humidity<0) = 0;

Absolute_Humidity = zeros(size(Model.T));
Model.WV = Absolute_Humidity;

%make water vapor the same dimentions as data
%===interpolate to new time grid===
%%%%abs_humid = interp2(double(time_Absolute_Humidity),double(range_Absolute_Humidity),double(Absolute_Humidity),Time.ts,Range.rm,'nearest',nan); % interpolate water vapor to range and time to match other data
%===integrate to new time grid===
%%%%[abs_humidintSum,abs_humidbins] = intSum(double(Absolute_Humidity),double(time_Absolute_Humidity),Counts.NBinsPython,Time.ts);
%shorten to new rm
%%%abs_humid = abs_humid(1:Range.i_range,:);
% % for ii = 1:Time.i_time
% %     abs_humid(:,ii) = interp1(double(range_Absolute_Humidity),abs_humidintSum(:,ii),Range.rm);
% % end
% % abs_humid = fillmissing(abs_humid,'nearest',1);
% % abs_humid = fillmissing(abs_humid,'nearest',2);
% % 
% % WV_intp = abs_humid / 1000 / Constant.mWV;        %[molecules/m^3] water vapor number density
% % 
% % % Fill in missing lower data with repeated values
% % WV_cut = 9;
% % WV_intp_fill = [NaN((WV_cut),Time.i_time); WV_intp(WV_cut+1:end,:)];
% % WV_intp_fill = fillmissing(WV_intp_fill,'nearest');
% % 
% % Model.WV = WV_intp_fill;                      %[molecules/m^3] water vapor number density
% % 
% % Model.WV = WV_intp;
% % %WV = zeros(size(WV));

%%
%=============================
%== ChristmanField sonde =========
%===========================

disp('Loading Sonde data')
[sonde_datetime,sondeStruc] =  ChristmanFieldradiosonde(pathSonde,span_days);
rm_sgp = cell(1,numel(sonde_datetime));
T_sonde_int = rm_sgp;
P_sonde_int = rm_sgp;
WV_sonde_int = rm_sgp;
rm_sonde_int = rm_sgp;
for i = 1:numel(sonde_datetime) % Loop over number of sondes in time period
    if isdatetime(sonde_datetime(i)) %== Check if sonde exists
        % ===Subtract first range value (site elevation) from whole vector
        %rm_sgp{i} = sondeStruc(i).Height - sondeStruc(i).Height(1);
        rm_sgp{i} = sondeStruc(i).Height - 1571.9;

        rm_sgp{i} = sondeStruc(i).Height - 1571.9;
        %===convert to same units====
        sondeStruc(i).P = sondeStruc(i).P./1013.25;%atm

        Sonde.Psurf(:,i) = sondeStruc(i).P(1);
        Sonde.Tsurf(:,i) = sondeStruc(i).T(1);
        % ==Collect radiosonde surface measurements==
%         T_sgp_surf(i) = sondeStruc(i).T(1);
%         P_sgp_surf(i) = sondeStruc(i).P(1);
        % ==Custom interpolation function==
        [T_sonde_int{i},P_sonde_int{i},WV_sonde_int{i},rm_sonde_int{i}] = interp_sonde2(sondeStruc(i).T,sondeStruc(i).P,sondeStruc(i).WV,rm_sgp{i},Range.rangeBin);  
        

        [rm_sgp{i},IA,~] = unique(rm_sgp{i});
        sondeStruc(i).T = sondeStruc(i).T(IA);
        sondeStruc(i).P = sondeStruc(i).P(IA);
        sondeStruc(i).WV = sondeStruc(i).WV(IA);
        T_sonde_int{i} = interp1(rm_sgp{i},sondeStruc(i).T,Range.rm,'makima',nan);
        P_sonde_int{i} = interp1(rm_sgp{i},sondeStruc(i).P,Range.rm,'makima',nan);
        WV_sonde_int{i} = interp1(rm_sgp{i},sondeStruc(i).WV,Range.rm,'makima',nan);
        
        if length(T_sonde_int{i})<Range.i_range % ==If sonde does not reach full lidar range
            disp('ran')
            Sonde.T_sonde(1:length(T_sonde_int{i}),i) = T_sonde_int{i};
            Sonde.T_sonde(length(T_sonde_int{i})+1:length(Range.rm),i)=nan(length(Range.rm)-length(T_sonde_int{i}),1);
            Sonde.P_sonde(1:length(T_sonde_int{i}),i) = P_sonde_int{i};
            Sonde.P_sonde(length(P_sonde_int{i})+1:length(Range.rm),i)=nan(length(Range.rm)-length(P_sonde_int{i}),1);
            Sonde.WV_sonde(1:length(T_sonde_int{i}),i) = WV_sonde_int{i};
            Sonde.WV_sonde(length(WV_sonde_int{i})+1:length(Range.rm),i)=nan(length(Range.rm)-length(WV_sonde_int{i}),1);
            Sonde.AbsHum(:,i) = Sonde.WV_sonde(:,i).*Constant.mWV*1000; %[g/m3]
        else
            Sonde.T_sonde(:,i) = T_sonde_int{i}(1:Range.i_range);
            Sonde.P_sonde(:,i) = P_sonde_int{i}(1:Range.i_range);
            Sonde.WV_sonde(:,i) = WV_sonde_int{i}(1:Range.i_range);
            Sonde.AbsHum(:,i) = Sonde.WV_sonde(:,i).*Constant.mWV*1000; %[g/m3]
        end
        %===interp sonde time
        [rm_sgp{i},IA,~] = unique(rm_sgp{i});
        sonde_time(1:length(rm_sonde_int{i}),i) = interp1(rm_sgp{i},sondeStruc(i).time(IA),rm_sonde_int{i})';  
        if length(sonde_time) < Range.i_range
            sonde_time = [sonde_time; sonde_time(end).*ones(Range.i_range-length(sonde_time),1)];
        end
         %Find index of sonde in time vector
         for j = 1:Range.i_range
            [~, Sonde.sonde_ind(j,i)]=min(abs(sonde_datetime(i)+seconds(sonde_time(j,i))-Time.date_ts));
         end
    else
        Sonde.sonde_ind = [];
        Sonde.T_sonde = nan(Range.i_range,1);
        Sonde.P_sonde = nan(Range.i_range,1);
        Sonde.WV_sonde = nan(Range.i_range,1);
        Sonde.AbsHum = nan(Range.i_range,1);
    end 
end


%=== set model WV to sonde interpolation
if ~isempty(Sonde.sonde_ind)
    if length(Sonde.sonde_ind(1,:))==1 %==if only one sonde
        Model.WV = ones(size(Model.WV)).*Sonde.WV_sonde;
        Model.WV = fillmissing(Model.WV,'linear',1);
        Model.WV = fillmissing(Model.WV,'linear',2);
    else
        Model.WV = interp2(Time.ts(Sonde.sonde_ind(1,:)),Range.rm,Sonde.WV_sonde,Time.ts,Range.rm,'linear');
        Model.WV = fillmissing(Model.WV,'linear',1);
        Model.WV = fillmissing(Model.WV,'linear',2);
    end
end

%%
%-----BSR
% %===integrate to new time grid===
% [BSRintSum,BSRbins] = intSum(double(Backscatter_Ratio),double(time_Backscatter_Ratio),Counts.NBinsPython,Time.ts);
% [Aerosol_Backscatter_CoefficientintSum,BSRbins] = intSum(double(Aerosol_Backscatter_Coefficient),double(time_Backscatter_Ratio),Counts.NBinsPython,Time.ts);
% for ii = 1:Time.i_time
%     HSRL.BSR(:,ii) = interp1(double(range_Backscatter_Ratio),BSRintSum(:,ii),Range.rm);
%     HSRL.Ba(:,ii) = interp1(double(range_Backscatter_Ratio),Aerosol_Backscatter_CoefficientintSum(:,ii),Range.rm);
% end
% HSRL.BSR = fillmissing(HSRL.BSR,'nearest',1);
% HSRL.BSR = fillmissing(HSRL.BSR,'nearest',2);
% 
% HSRL.Ba = fillmissing(HSRL.Ba,'nearest',1);
% HSRL.Ba = fillmissing(HSRL.Ba,'nearest',2);
% 
% HSRL.Bm = HSRL.Ba./(HSRL.BSR-1);%'m^(-1) sr^(-1)'


%%
% --- O2 Background Subtraction ---
if strcmp(Options.MPDname,'05') || strcmp(Options.MPDname,'03')
    Counts.bg_o2on = round(mean(Data.MCS.Channel9.Data(end-20:end,:)));% Take mean of last data points
    Counts.o2on_bgsub = Data.MCS.Channel9.Data - Counts.bg_o2on;       % Background subtracted
    %Counts.o2on_bgsub(Counts.o2on_bgsub < 0) = 0;         % Minimum of zero
    
    Counts.bg_o2off = round(mean(Data.MCS.Channel1.Data(end-20:end,:)));% Take mean of last data points
    Counts.o2off_bgsub = Data.MCS.Channel1.Data - Counts.bg_o2off;       % Background subtracted
    %Counts.o2off_bgsub(Counts.o2off_bgsub < 0) = 0;         % Minimum of zero
    
    Counts.bg_o2on_mol = round(mean(Data.MCS.Channel10.Data(end-20:end,:)));% Take mean of last data points
    Counts.o2on_bgsub_mol = Data.MCS.Channel10.Data - Counts.bg_o2on_mol;       % Background subtracted
   % Counts.o2on_bgsub_mol(Counts.o2on_bgsub_mol < 0) = 0;         % Minimum of zero
    
    Counts.bg_o2off_mol = round(mean(Data.MCS.Channel2.Data(end-20:end,:)));% Take mean of last data points
    Counts.o2off_bgsub_mol = Data.MCS.Channel2.Data - Counts.bg_o2off_mol;       % Background subtracted
    %Counts.o2off_bgsub_mol(Counts.o2off_bgsub_mol < 0) = 0;         % Minimum of zero
    
    Counts.bg_wvon = round(mean(Data.MCS.Channel8.Data(end-20:end,:)));% Take mean of last data points
    Counts.wvon_bgsub = Data.MCS.Channel8.Data - Counts.bg_wvon;       % Background subtracted
    %Counts.wvon_bgsub(Counts.wvon_bgsub < 0) = 0;         % Minimum of zero
    
    Counts.bg_wvoff = round(mean(Data.MCS.Channel0.Data(end-20:end,:)));% Take mean of last data points
    Counts.wvoff_bgsub = Data.MCS.Channel0.Data - Counts.bg_wvoff;       % Background subtracted
    %Counts.wvoff_bgsub(Counts.wvoff_bgsub < 0) = 0;         % Minimum of zero

elseif strcmp(Options.MPDname,'01')

    Counts.bg_o2off = mean(Data.MCS.Channel9.Data(end-20:end,:));% Take mean of last data points
    Counts.o2off_bgsub = Data.MCS.Channel9.Data - Counts.bg_o2off;       % Background subtracted
    Counts.o2off_bgsub(Counts.o2off_bgsub < 0) = 0;         % Minimum of zero
    
    Counts.bg_o2on = mean(Data.MCS.Channel1.Data(end-20:end,:));% Take mean of last data points
    Counts.o2on_bgsub = Data.MCS.Channel1.Data - Counts.bg_o2on;       % Background subtracted
    Counts.o2on_bgsub(Counts.o2on_bgsub < 0) = 0;         % Minimum of zero
    
    Counts.bg_o2off_mol = mean(Data.MCS.Channel10.Data(end-20:end,:));% Take mean of last data points
    Counts.o2off_bgsub_mol = Data.MCS.Channel10.Data - Counts.bg_o2off_mol;       % Background subtracted
    Counts.o2off_bgsub_mol(Counts.o2off_bgsub_mol < 0) = 0;         % Minimum of zero
    
    Counts.bg_o2on_mol = mean(Data.MCS.Channel2.Data(end-20:end,:));% Take mean of last data points
    Counts.o2on_bgsub_mol = Data.MCS.Channel2.Data - Counts.bg_o2on_mol;       % Background subtracted
    Counts.o2on_bgsub_mol(Counts.o2on_bgsub_mol < 0) = 0;         % Minimum of zero
    
    Counts.bg_wvon = mean(Data.MCS.Channel8.Data(end-20:end,:));% Take mean of last data points
    Counts.wvon_bgsub = Data.MCS.Channel8.Data - Counts.bg_wvon;       % Background subtracted
    Counts.wvon_bgsub(Counts.wvon_bgsub < 0) = 0;         % Minimum of zero
    
    Counts.bg_wvoff = mean(Data.MCS.Channel0.Data(end-20:end,:));% Take mean of last data points
    Counts.wvoff_bgsub = Data.MCS.Channel0.Data - Counts.bg_wvoff;       % Background subtracted
    Counts.wvoff_bgsub(Counts.wvoff_bgsub < 0) = 0;         % Minimum of zero
end

%%
%integrate to new range
inc = 1;
for ii = 1:2:length(Range.rm_raw_o2)
    o2on_intp2(inc,:) = (Counts.o2on_bgsub(ii,:)+Counts.o2on_bgsub(ii+1,:));
    o2off_intp2(inc,:) = (Counts.o2off_bgsub(ii,:)+Counts.o2off_bgsub(ii+1,:));
    o2on_intp2_mol(inc,:) = (Counts.o2on_bgsub_mol(ii,:)+Counts.o2on_bgsub_mol(ii+1,:));
    o2off_intp2_mol(inc,:) = (Counts.o2off_bgsub_mol(ii,:)+Counts.o2off_bgsub_mol(ii+1,:));
    wvon_intp2(inc,:) = (Counts.wvon_bgsub(ii,:)+Counts.wvon_bgsub(ii+1,:));
    wvoff_intp2(inc,:) = (Counts.wvoff_bgsub(ii,:)+Counts.wvoff_bgsub(ii+1,:));
    inc = inc+1;
end
Counts.o2on_bgsub = o2on_intp2;
Counts.o2off_bgsub = o2off_intp2;
Counts.o2on_bgsub_mol = o2on_intp2_mol;
Counts.o2off_bgsub_mol = o2off_intp2_mol;
Counts.wvon_bgsub = wvon_intp2;
Counts.wvoff_bgsub = wvoff_intp2;

Range.rm_raw_o2 = Range.rm_raw_o2(1:2:end)+Range.rangeBin./2;

%%

% Interpolating to shorter range vector
% Counts.o2on_noise = interp2(Time.ts,Range.rm_raw_o2,Counts.o2on_bgsub,Time.ts,Range.rm);
% Counts.o2on_noise_mol = interp2(Time.ts,Range.rm_raw_o2,Counts.o2on_bgsub_mol,Time.ts,Range.rm);
% Counts.o2off_noise = interp2(Time.ts,Range.rm_raw_o2,Counts.o2off_bgsub,Time.ts,Range.rm);
% Counts.o2off_noise_mol = interp2(Time.ts,Range.rm_raw_o2,Counts.o2off_bgsub_mol,Time.ts,Range.rm);
% Counts.wvon_noise = interp2(Time.ts,Range.rm_raw_o2,Counts.wvon_bgsub,Time.ts,Range.rm);
% Counts.wvoff_noise = interp2(Time.ts,Range.rm_raw_o2,Counts.wvoff_bgsub,Time.ts,Range.rm);


Counts.o2on_noise = Counts.o2on_bgsub(1:80,:);
Counts.o2on_noise_mol = Counts.o2on_bgsub_mol(1:80,:);
Counts.o2off_noise = Counts.o2off_bgsub(1:80,:);
Counts.o2off_noise_mol = Counts.o2off_bgsub_mol(1:80,:);
Counts.wvon_noise = Counts.wvon_bgsub(1:80,:);
Counts.wvoff_noise = Counts.wvoff_bgsub(1:80,:);

% Counts.o2on_noise = Counts.o2on_bgsub(7:86,:);
% Counts.o2on_noise_mol = Counts.o2on_bgsub_mol(7:86,:);
% Counts.o2off_noise = Counts.o2off_bgsub(7:86,:);
% Counts.o2off_noise_mol = Counts.o2off_bgsub_mol(7:86,:);
% Counts.wvon_noise = Counts.wvon_bgsub(7:86,:);
% Counts.wvoff_noise = Counts.wvoff_bgsub(7:86,:);

Counts.NBins = Data.MCS.Channel0.NBins.*2;

%%
%-----Create Spectrum vectors----

%lambda_online = interp1(Options.TimeGrid,Data.Laser.O2Online.WavelengthActual,Time.ts/60/60);
%lambda_offline = interp1(Options.TimeGrid,Data.Laser.O2Offline.WavelengthActual,Time.ts/60/60);

Spectrum.WavemeterOffset = 0.00024;%[nm]
%Spectrum.WavemeterOffset = 0;%[nm]


Spectrum.lambda_online = 769.7958 *ones(size(Time.ts));
Spectrum.lambda_offline = 770.1085 *ones(size(Time.ts));

% Spectrum.lambda_wvon = 828.1959 *ones(size(Time.ts));
% Spectrum.lambda_wvoff = 828.2951 *ones(size(Time.ts));
% Spectrum.lambda_wvon = 828.187 *ones(size(Time.ts));
% Spectrum.lambda_wvoff = 828.2951 *ones(size(Time.ts));
Spectrum.lambda_wvon = fillmissing(filloutliers(Data.Laser.WVOnline.WavelengthActual,'linear','movmedian',5),'linear');
Spectrum.lambda_wvoff = fillmissing(filloutliers(Data.Laser.WVOffline.WavelengthActual,'linear','movmedian',5),'linear');

Spectrum.lambda_wvon = double(Spectrum.lambda_wvon); %single values mess up conversion to wavenumber
Spectrum.lambda_wvoff = double(Spectrum.lambda_wvoff);

wvlambdaCentralon = 828.1965;
wvnuCentralon = 10^7/wvlambdaCentralon;
%wvlambdaCentraloff = 828.2957;

Spectrum.nu_online = 10^7./Spectrum.lambda_online;                    %[cm-1] Online wavenumber
Spectrum.nu_offline = 10^7./Spectrum.lambda_offline;                  %[cm-1] Offline wavenumber

Spectrum.nu_wvon = 10^7./Spectrum.lambda_wvon;                    %[cm-1] Online wavenumber
Spectrum.nu_wvoff = 10^7./Spectrum.lambda_wvoff;                  %[cm-1] Offline wavenumber

nuMin = Spectrum.nu_online-0.334;                                 %[cm-1] Scan lower bound
nuMax = Spectrum.nu_online+0.334;                                 %[cm-1] Scan upper bound
Spectrum.nuBin = 0.00222;                                    %[cm-1] Scan increment
nu_scan = (nuMin:Spectrum.nuBin:nuMax);                      %[cm-1](1 x nu) Scan vector

nuwvMin = wvnuCentralon-0.334;                                 %[cm-1] Scan lower bound
nuwvMax = wvnuCentralon+0.334;                                 %[cm-1] Scan upper bound
Spectrum.nuBin = 0.00222;                                    %[cm-1] Scan increment
nu_scanwv = (nuwvMin:Spectrum.nuBin:nuwvMax);                      %[cm-1](1 x nu) Scan vector

nuMin_off = Spectrum.nu_offline-0.334;                                 %[cm-1] Scan lower bound
nuMax_off = Spectrum.nu_offline+0.334;                                 %[cm-1] Scan upper bound
nu_scan_off = (nuMin_off:Spectrum.nuBin:nuMax_off);

Spectrum.nu_scan_3D_short = permute(nu_scan, [3 1 2]);       %[cm-1] putting scan in third dimension
Spectrum.nu_scan_3D_short_off = permute(nu_scan_off, [3 1 2]);       %[cm-1] putting scan in third dimension

Spectrum.nu_scanwv_3D_short = permute(nu_scanwv, [3 1 2]);       %[cm-1] putting scan in third dimension

Spectrum.lambda_scan_3D_short = 10^7./Spectrum.nu_scan_3D_short;
Spectrum.lambda_scan_3D_short_off = 10^7./Spectrum.nu_scan_3D_short_off;
Spectrum.i_scan_3D_short = length(Spectrum.nu_scan_3D_short);         %[none] length of scan vector

Spectrum.lambda_scanwv_3D_short = 10^7./Spectrum.nu_scan_3D_short;

%add wavemeter offset
Spectrum.lambda_online = 769.7958 *ones(size(Time.ts))+Spectrum.WavemeterOffset;
Spectrum.lambda_offline = 770.1085 *ones(size(Time.ts))+Spectrum.WavemeterOffset;
Spectrum.lambda_wvon = fillmissing(filloutliers(Data.Laser.WVOnline.WavelengthActual,'linear','movmedian',5),'linear');
Spectrum.lambda_wvoff = fillmissing(filloutliers(Data.Laser.WVOffline.WavelengthActual,'linear','movmedian',5),'linear');
Spectrum.lambda_wvon = double(Spectrum.lambda_wvon)+Spectrum.WavemeterOffset; %single values mess up conversion to wavenumber
Spectrum.lambda_wvoff = double(Spectrum.lambda_wvoff)+Spectrum.WavemeterOffset;
Spectrum.nu_online = 10^7./Spectrum.lambda_online;                    %[cm-1] Online wavenumber
Spectrum.nu_offline = 10^7./Spectrum.lambda_offline;                  %[cm-1] Offline wavenumber
Spectrum.nu_wvon = 10^7./Spectrum.lambda_wvon;                    %[cm-1] Online wavenumber
Spectrum.nu_wvoff = 10^7./Spectrum.lambda_wvoff;                  %[cm-1] Offline wavenumber


[~,Spectrum.online_index] = min(abs(Spectrum.nu_online - Spectrum.nu_scan_3D_short),[],3);%finding index of online wavenumber
[~,Spectrum.offline_index] = min(abs(Spectrum.nu_offline - Spectrum.nu_scan_3D_short_off),[],3);%finding index of online wavenumber

[~,Spectrum.online_indexwv] = min(abs(Spectrum.nu_wvon - Spectrum.nu_scanwv_3D_short),[],3);%finding index of online wavenumber

Model.absorption = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_online,Model.WV); %[m-1] Funcrtion to calculate theoretical absorption
%Model.transmission = exp(-cumtrapz(Range.rm,Model.absorption));
Model.absorption_off = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_offline,Model.WV); %[m-1] Funcrtion to calculate theoretical absorption

%%
%Count filtering
k = ones(Options.oversample,Options.t_avg)./(Options.oversample*Options.t_avg);     % Kernel

Counts.o2on = filter2(k,Counts.o2on_noise,'same');
Counts.o2on_mol = filter2(k,Counts.o2on_noise_mol,'same');
Counts.o2off = filter2(k,Counts.o2off_noise,'same');
Counts.o2off_mol = filter2(k,Counts.o2off_noise_mol,'same');

Counts.wvon = filter2(k,Counts.wvon_noise,'same');
Counts.wvoff = filter2(k,Counts.wvoff_noise,'same');

Counts.o2on = Counts.o2on_noise;
Counts.o2on_mol = Counts.o2on_noise_mol;
Counts.o2off = Counts.o2off_noise;
Counts.o2off_mol = Counts.o2off_noise_mol;

Counts.wvon = Counts.wvon_noise;
Counts.wvoff = Counts.wvoff_noise;

%%
% % load('overlapCorrMSU6_18_21.mat','overlapCorr')
% load('overlapCorrMSU6_21_21.mat','overlapCorr')
% %Counts.o2off = Counts.o2off.*overlapCorr;

% load('overlapCorrMSU6_24_21.mat','overlapCorr')
% load('overlapCorrMSU6_24_21_1000.mat','overlapCorr')
%  [~,constAlt]=min(abs(2000-Range.rm));
%  Counts.o2off(1:constAlt,:) = Counts.o2off(1:constAlt,:).*overlapCorr(1:constAlt,:);

% load('BoulderlnOonOoff6_16to6_29_21.mat','lnOonOoff')
% lnOonOoff = smoothdata(lnOonOoff,2,'movmean',10);
%   [~,constAlt]=min(abs(500-Range.rm));
%   Counts.o2off(1:constAlt,:) = Counts.o2off(1:constAlt,:).*lnOonOoff(1:constAlt,:);
%%

%====== Calucate any appy optimal filtering based on Poisson thinning ====
tic
%%%Counts = poissonThin(Counts);
toc
%%

% Pulse Decon
% pulse = [0.5;0.5];
% onConv = conv2(Counts.o2on,pulse);
% Counts.o2on = 2*Counts.o2on - onConv(2:end,:);
% Counts.o2on(Counts.o2on<0)=0;
% offConv = conv2(Counts.o2off,pulse);
% Counts.o2off = 2*Counts.o2off - offConv(2:end,:);
% Counts.o2off(Counts.o2off<0)=0;
% 
% on_molConv = conv2(Counts.o2on_mol,pulse);
% Counts.o2on_mol = 2*Counts.o2on_mol - on_molConv(2:end,:);
% Counts.o2on_mol(Counts.o2on_mol<0)=0;
% off_molConv = conv2(Counts.o2off_mol,pulse);
% Counts.o2off_mol = 2*Counts.o2off_mol - off_molConv(2:end,:);
% Counts.o2off_mol(Counts.o2off_mol<0)=0;

%%
%===== Calucation Model absorption for radiosondes =======
for i=1:numel(sonde_datetime) 
        if isdatetime(sonde_datetime(i)) %Check if there are any sondes
           % Sonde.absorption_sonde{i} = diag(absorption_O2_770_model(Sonde.T_sonde(:,i),Sonde.P_sonde(:,i),Spectrum.nu_online(Sonde.sonde_ind(:,i)),Model.WV(:,Sonde.sonde_ind(:,i)))); %[m-1] Funcrtion to calculate theoretical absorption
           Sonde.absorption_sonde{i} = absorption_O2_770_model(Sonde.T_sonde(:,i),Sonde.P_sonde(:,i),Spectrum.nu_online(Sonde.sonde_ind(1,i)),Sonde.WV_sonde(:,i)); %[m-1] Funcrtion to calculate theoretical absorption
            Sonde.trasmission_sonde{i} = exp(-cumtrapz(Range.rm,Sonde.absorption_sonde{i})); %O2 transmission
        else
            Sonde.absorption_sonde{i} = nan(Range.i_range,1);
            Sonde.trasmission_sonde{i} = nan(Range.i_range,1);
            Sonde.T_sonde = nan(Range.i_range,1);
            Sonde.P_sonde = nan(Range.i_range,1);
        end
end
%%
Data.Thermocouple.TSOA.Temperature = nan(size(Data.Thermocouple.InsideCell.Temperature));
Data.Thermocouple.OutsideCell.Temperature = nan(size(Data.Thermocouple.InsideCell.Temperature));

%%
disp('Calculating HSRL')


if strcmp(Options.MPDname,'03')
            Counts.Nc_on = Counts.o2on;%.*Counts.NBins;
    Counts.Nc_off = Counts.o2off;%.*Counts.NBins;
    Counts.Nm_on = Counts.o2on_mol;%.*Counts.NBins;
    Counts.Nm_off = Counts.o2off_mol;%.*Counts.NBins;
    [HSRL] = backscatterRetrievalMPD03(Counts, Model, Spectrum);

else

    LidarData.Range = Range.rm;
    LidarData.Time = Time.ts;
    LidarData.OfflineCombinedTotalCounts = Counts.o2off;
    LidarData.OfflineMolecularTotalCounts = Counts.o2off_mol;
    WeatherData.Temperature = Model.T;
    WeatherData.Pressure = Model.P;
    if span_days(1) <  datetime(2021,12,1,1,0,0,'TimeZone','UTC')
        [LidarData]=BackscatterRetrievalBoulder062021(LidarData,WeatherData);
    else
        [LidarData]=BackscatterRetrievalBoulder031122(LidarData,WeatherData);
    end
    HSRL.BSR = LidarData.UnmaskedBackscatterRatio;
    HSRL.Ba = LidarData.UnmaskedAerosolBackscatterCoefficient;
    HSRL.Bm = LidarData.MolecularBackscatterCoefficient;
    %Ba828 = LidarData.UnmaskedAerosolBackscatterCoefficient828;
    %Bm828 = LidarData.MolecularBackscatterCoefficient828;
    HSRL.BSR828 =LidarData.UnmaskedBackscatterRatio828;
end




% 
%     LidarData.OfflineCombinedTotalCounts = Counts.foff-Counts.foff_bg;
%     LidarData.OfflineMolecularTotalCounts = Counts.foff_mol-Counts.foff_mol_bg;
%     [LidarData]=BackscatterRetrievalBoulder062021(LidarData,WeatherData);
%     HSRL.BSRf = LidarData.UnmaskedBackscatterRatio;
% 
%     LidarData.OfflineCombinedTotalCounts = Counts.goff-Counts.foff_bg;
%     LidarData.OfflineMolecularTotalCounts = Counts.goff_mol-Counts.foff_mol_bg;
%     [LidarData]=BackscatterRetrievalBoulder062021(LidarData,WeatherData);
%     HSRL.BSRg = LidarData.UnmaskedBackscatterRatio;


    %HSRL.BSR = interp2(time_Backscatter_Ratio,range_Backscatter_Ratio,Backscatter_Ratio,Time.ts,Range.rm);
%%
%[Model] = modelCounts(Counts,Model,HSRL,Time,Range,Spectrum);