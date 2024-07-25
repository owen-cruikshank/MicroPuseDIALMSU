%load MSU data files
function [DataStructure2, Options] = loadMSUNETcdf(spanDays,Options)
%File: loadMSUNETcdf.m
%Date: 05/21/2020
%Author: Owen Cruikshank
%Inputs:
%   -spanDays:(datetime) Day to read data
%   -path:(string) string of file path were RSync folder data is stored
%
%Outputs:
%   -DataStructure2: (structure) All data read in
%   -Options: (structure) Data options
%%
%Loop over multiple days
for i=1:length(spanDays)

    % Convert datetime vector to YYYYMMDD sting
    [y,m,d]=ymd(spanDays(i));
    y=num2str(y);
    m=num2str(m);
    d=num2str(d);
    if length(m)<2 
        m=['0' m];
    end
    if length(d)<2 
        d=['0' d];
    end
    date = [y m d];

    %% Defining options
    Options.Date{i}          = date;
end
    %Options.System        = 'DIAL05';
    Options.InterpMethod  = 'linear';
    Options.Extrapolation = 'extrap';

    % Time grid to average bins over
    % To read in data starting after a 00:00:00 the starting time must be
    % changed
    Options.TimeGrid      = (0.1:Options.intTime:60*(24*length(spanDays)-1/60/60))./60;
    Options.TimeGrid      = (0.001:Options.intTime:60*(24*length(spanDays)-1/60/60))./60;

    %% Loading filepaths
    Paths.Code = pwd;
    
    for i=1:length(spanDays)
        Paths.Data{:,i} = fullfile(Options.path,Options.Date{i});
    end

    %% Defining data to read
    if strcmp(Options.MPDname,'05') || strcmp(Options.MPDname,'01') || strcmp(Options.MPDname,'03')
    DataTypes = {'Etalon*.nc';'LL*.nc';'MCS*.nc';'Power*.nc'; 'HKeep*.nc'; 'UPS*.nc'; 'WS*.nc'};
    DataNames = {'Etalon';'Laser';'MCS';'Power';'Thermocouple';'UPS';'WS'};
    FileVarNames  = {{'time';'Temperature';'TempDiff';'IsLocked';'EtalonNum'};
                     {'time';'Wavelength';'WaveDiff';'IsLocked';'TempDesired';'TempMeas';'Current';'SeedPower';'LaserName'};
                     {'time';'ProfilesPerHist';'Channel';'nsPerBin';'NBins';'Data';'RTime'};
                     {'time';'RTime';'Power';'AccumEx';'Demux'};
                     {'time';'Temperature'};
                     {'time';'BatteryNominal';'BatteryReplace';'BatteryInUse';'BatteryLow';'BatteryCapacity';'BatteryTimeLeft';'UPSTemperature';'HoursOnBattery'};
                     {'time';'Temperature';'RelHum';'Pressure';'AbsHum'}};  
    CodeVarNames  = {{'TimeStamp';'TemperatureActual';'-';'-';'Type'};
                     {'TimeStamp';'WavelengthActual';'-';'Locked';'TemperatureDesired';'TemperatureActual';'-';'-';'Type'};
                     {'TimeStamp';'ProfilesPerHistogram';'ChannelNum';'RangeResolution';'-';'-';'-'};
                     {'TimeStamp';'-';'LaserPower';'-';'-'};
                     {'TimeStamp';'Temperature'};
                     {'TimeStamp';'BatteryNominal';'BatteryReplace';'BatteryInUse';'BatteryLow';'BatteryCapacity';'BatteryTimeLeft';'UPSTemperature';'HoursOnBattery'};
                     {'TimeStamp';'Temperature';'RelHum';'Pressure';'AbsHum'}}; 
%     VarDesiredTypes = {{'-';'-';'-';'-';'String'};
%                        {'-';'-';'-';'-';'-';'-';'-';'-';'String'};
%                        {'-';'-';'-';'-';'-';'-';'-'};
%                        {'-';'-';'-';'-';'-'};
%                        {'-';'-'};
%                        {'-';'-';'-';'-';'-';'-';'-';'-';'-'};
%                        {'-';'-';'-';'-';'-'}};

        VarDesiredTypes = {{'-';'-';'-';'-';'String'};
                       {'-';'Double';'-';'-';'-';'-';'-';'-';'String'};
                       {'Double';'Double';'Double';'Double';'Double';'Double';'Double'};
                       {'-';'-';'-';'-';'-'};
                       {'-';'-'};
                       {'-';'-';'-';'-';'-';'-';'-';'-';'-'};
                       {'Double';'Double';'Double';'Double';'Double'}};
    elseif strcmp(Options.MPDname,'000')
            DataTypes = {'Etalon*.nc';'LL*.nc';'MCS*.nc';'Power*.nc'; 'UPS*.nc'};
    DataNames = {'Etalon';'Laser';'MCS';'Power';'UPS'};
    FileVarNames  = {{'time';'Temperature';'TempDiff';'IsLocked';'EtalonNum'};
                     {'time';'Wavelength';'WaveDiff';'IsLocked';'TempDesired';'TempMeas';'Current';'SeedPower';'LaserName'};
                     {'time';'ProfilesPerHist';'Channel';'nsPerBin';'NBins';'Data';'RTime'};
                     {'time';'RTime';'Power';'AccumEx';'Demux'};
                     {'time';'BatteryNominal';'BatteryReplace';'BatteryInUse';'BatteryLow';'BatteryCapacity';'BatteryTimeLeft';'UPSTemperature';'HoursOnBattery'}};  
    CodeVarNames  = {{'TimeStamp';'TemperatureActual';'-';'-';'Type'};
                     {'TimeStamp';'WavelengthActual';'-';'Locked';'TemperatureDesired';'TemperatureActual';'-';'-';'Type'};
                     {'TimeStamp';'ProfilesPerHistogram';'ChannelNum';'RangeResolution';'-';'-';'-'};
                     {'TimeStamp';'-';'LaserPower';'-';'-'};
                     {'TimeStamp';'BatteryNominal';'BatteryReplace';'BatteryInUse';'BatteryLow';'BatteryCapacity';'BatteryTimeLeft';'UPSTemperature';'HoursOnBattery'}}; 
    VarDesiredTypes = {{'-';'-';'-';'-';'String'};
                       {'-';'-';'-';'-';'-';'-';'-';'-';'String'};
                       {'-';'-';'-';'-';'-';'-';'-'};
                       {'-';'-';'-';'-';'-'};
                       {'-';'-';'-';'-';'-';'-';'-';'-';'-'}};
    else
    DataTypes = {'Etalon*.nc';'LL*.nc';'MCS*.nc';'Power*.nc'; 'HKeep*.nc'; 'UPS*.nc'};
    DataNames = {'Etalon';'Laser';'MCS';'Power';'Thermocouple';'UPS'};
    FileVarNames  = {{'time';'Temperature';'TempDiff';'IsLocked';'EtalonNum'};
                     {'time';'Wavelength';'WaveDiff';'IsLocked';'TempDesired';'TempMeas';'Current';'SeedPower';'LaserName'};
                     {'time';'ProfilesPerHist';'Channel';'nsPerBin';'NBins';'Data';'RTime'};
                     {'time';'RTime';'Power';'AccumEx';'Demux'};
                     {'time';'Temperature'};
                     {'time';'BatteryNominal';'BatteryReplace';'BatteryInUse';'BatteryLow';'BatteryCapacity';'BatteryTimeLeft';'UPSTemperature';'HoursOnBattery'}};  
    CodeVarNames  = {{'TimeStamp';'TemperatureActual';'-';'-';'Type'};
                     {'TimeStamp';'WavelengthActual';'-';'Locked';'TemperatureDesired';'TemperatureActual';'-';'-';'Type'};
                     {'TimeStamp';'ProfilesPerHistogram';'ChannelNum';'RangeResolution';'-';'-';'-'};
                     {'TimeStamp';'-';'LaserPower';'-';'-'};
                     {'TimeStamp';'Temperature'};
                     {'TimeStamp';'BatteryNominal';'BatteryReplace';'BatteryInUse';'BatteryLow';'BatteryCapacity';'BatteryTimeLeft';'UPSTemperature';'HoursOnBattery'}}; 
    VarDesiredTypes = {{'-';'-';'-';'-';'String'};
                       {'-';'-';'-';'-';'-';'-';'-';'-';'String'};
                       {'-';'-';'-';'-';'-';'-';'-'};
                       {'-';'-';'-';'-';'-'};
                       {'-';'-'};
                       {'-';'-';'-';'-';'-';'-';'-';'-';'-'}};
    end


    % Populating cells with "-" variables with their parner
    for m=1:1:size(CodeVarNames,1)
       for n=1:1:size(CodeVarNames{m,1},1) 
          if strcmp(CodeVarNames{m,1}{n},'-')
              CodeVarNames{m,1}{n} = FileVarNames{m,1}{n};
          end
          if strcmp(VarDesiredTypes{m,1}{n},'-')
              %%VarDesiredTypes{m,1}{n} = 'Double';
              VarDesiredTypes{m,1}{n} = 'Native';
          end
       end
    end   
    % Making a single cell array to pass to the loading function
    for m=1:1:size(DataTypes,1)
        ToLoad{m,1} = DataTypes{m,1}; %#ok<*SAGROW>
        ToLoad{m,2} = FileVarNames{m,1};
        ToLoad{m,3} = CodeVarNames{m,1};
        ToLoad{m,4} = VarDesiredTypes{m,1};
    end     

    clear m n CodeVarNames DataTypes FileVarNames  VarDesiredTypes Data2 DataStructure

    %% Loading data
    for i = 1:length(spanDays)
        disp(spanDays(i))
        RawData = ReadMPDData(ToLoad,Paths.Code,Paths.Data{:,i});

        for m=1:size(RawData,1)
            for p=1:size(RawData{m,1},1)
                if isequal(ToLoad{m,3}{p,1},'TimeStamp')
                        RawData{m,1}{p,1} = RawData{m,1}{p,1} + 24*(i-1);
                end
            end
        end

    %% Parsing data
    disp('Parsing Data')
    if strcmp(Options.MPDname,'05') || strcmp(Options.MPDname,'01') || strcmp(Options.MPDname,'03')
        DemuxNames = {{'O2Etalon';'WVEtalon'};
             {'O2Online';'O2Offline';'WVOnline';'WVOffline'};
             {'Channel0';'Channel1';'Channel2';'Channel8';'Channel9';'Channel10'};
             {'O2Online';'O2Offline'};
             {'RoomTemp','HVACReturn','HVACSource','Window','WindowHeaterFan','VWEtalonHeatSink','HSRLEtalonHeatSink','InsideCell'}
             {'all'};
             {'all'}}; 
    % Channel labels
       Demux = {{'O2Etalon';'WVEtalon'};
             {'O2Online';'O2Offline';'WVOnline';'WVOffline'};
             {0;1;2;8;9;10};
             {'O2Online';'Unknown';'Unknown';'Unknown';'Unknown';'Unknown';'O2Offline';'Unknown';'Unknown';'Unknown';'Unknown';'Unknown'};
             {1;2;3;4;5;6;7;8};
             {'UPS'};
             {'WS'}};
    % Column where names are stored in raw data
    DemuxCol = [5 9 3 5 0 0 0];
    
    
    elseif strcmp(Options.MPDname,'000')
             DemuxNames = {{'HSRLEtalon'};
             {'HSRL'};
             {'Channel0';'Channel2';'Channel8';'Channel10'};
             {'O2Online';'O2Offline'};
             {'all'};
             {'all'}}; 
      % Channel labels
       Demux = {{'HSRLEtalon'};
             {'HSRL'};
             {0;2;8;10};
             {'O2Online';'Unknown';'Unknown';'Unknown';'Unknown';'Unknown';'O2Offline';'Unknown';'Unknown';'Unknown';'Unknown';'Unknown'};
             {'UPS'};
             {'WS'}};
        % Column where names are stored in raw data
        DemuxCol = [5 9 3 5 0 0 0];
    else
               % Data names
        DemuxNames = {{'O2Etalon'};
                 {'O2Online';'O2Offline';'TWSOA'};
                 {'Channel0';'Channel2';'Channel8';'Channel10'};
                 {'O2Online';'O2Offline'};
                 {'InsideCell','OutsideCell','TSOA','RoomTemp'};
                 {'all'}}; 
    
        Demux = {{'O2Etalon'};
                 {'O2Online';'O2Offline';'Unknown'};
                 {0;2;8;10};
                 {'O2Online';'Unknown';'Unknown';'Unknown';'Unknown';'Unknown';'O2Offline';'Unknown';'Unknown';'Unknown';'Unknown';'Unknown'};
                 {1;2;3;4};
                 {'UPS'}};
                 
            % Column where names are stored in raw data
        DemuxCol = [5 9 3 5 0 0];
    end


    clear parsedData
    for m=1:size(RawData,1)             %Loop over variable names
        for n=1:size(Demux{m,1},1)      %Loop over demux options
           for p=1:size(RawData{m,1},1) %Loop over variable names
               % Power and thermocouples have different demuxing
               if strcmp(DataNames{m,1},'Power')
                   %%% Do nothing for now because we have no power measurements
                   if (p==1||p==2) && (n==1||n==7)% Check if time vector
                       if n==7
                           parsedData{m,1}{2,1}{p,1}(:,1) = RawData{m,1}{p,1}(:,1);
                       else
                           parsedData{m,1}{n,1}{p,1}(:,1) = RawData{m,1}{p,1}(:,1);
                       end
                   elseif n==1
                       parsedData{m,1}{n,1}{p,1}(:,1) = RawData{m,1}{p,1}(:,n);
                   elseif n==7
                       parsedData{m,1}{2,1}{p,1}(:,1) = RawData{m,1}{p,1}(:,n); 
                   end
               elseif strcmp(DataNames{m,1},'Thermocouple')
                   if p==1 % Check if time vector
                       parsedData{m,1}{n,1}{p,1}(:,1) = RawData{m,1}{p,1}(:,1);
                   else
                       parsedData{m,1}{n,1}{p,1}(:,1) = RawData{m,1}{p,1}(:,n); 
                   end
               elseif strcmp(DataNames{m,1},'UPS')
                   if p==1 % Check if time vector
                       parsedData{m,1}{n,1}{p,1}(:,1) = RawData{m,1}{p,1}(1:end,1);
                   else
                       parsedData{m,1}{n,1}{p,1}(:,1) = RawData{m,1}{p,1}(1:end,n); 
                   end
               elseif strcmp(DataNames{m,1},'WS')
                   if p==1 % Check if time vector
                       parsedData{m,1}{n,1}{p,1}(:,1) = RawData{m,1}{p,1}(1:end,1);
                   else
                       parsedData{m,1}{n,1}{p,1}(:,1) = RawData{m,1}{p,1}(1:end,n); 
                   end
               elseif strcmp(DataNames{m,1},'MCS')
                   % Create vector of T/F values where names equal data
                   demuxTrue = Demux{m,1}{n,1} == RawData{m,1}{DemuxCol(m),1}(:);
                   % Assign data corresponding to demux names to new cells
                   parsedData{m,1}{n,1}{p,1}(:,:) = RawData{m,1}{p,1}(demuxTrue,:);
               else % Demux for all other variables
                   demuxTrue = strcmp(Demux{m,1}{n,1},RawData{m,1}{DemuxCol(m),1}(:));
                   parsedData{m,1}{n,1}{p,1}(:,:) = RawData{m,1}{p,1}(demuxTrue,:);
               end
           end
        end
    end
    %%
    clear m n q p increment
    
    %% Pushing all 1d data to a constant grid
    disp('Interpolating data to grid')

TimeGrid = Options.TimeGrid(Options.TimeGrid>=((i-1)*24) & Options.TimeGrid<(i*24));

    for m=1:size(parsedData,1)        %Loop over variable names
        for n=1:size(parsedData{m,1},1) %Loop over demux options
            % Converting the cell array data to a structure
            DataStructure{m,1}{n,1} = cell2struct(parsedData{m,1}{n,1},ToLoad{m,3});
            % Interpolating 1d data or passing through MCS 2d data
            if strcmp(DataNames{m,1},'MCS') || strcmp(DataNames{m,1},'Power')
               % Checking if power time stamps are monotonically increasing
               A = find(diff(DataStructure{m,1}{n,1}.TimeStamp)<0) + 1;
               for p=1:1:length(A)
                    DataStructure{m,1}{n,1}.TimeStamp(A(p)) = DataStructure{m,1}{n,1}.TimeStamp(A(p)) +DataStructure{m,1}{n,1}.TimeStamp(A(p)-1) ;
               end
               % Average data over time grid
               % Initialize variables
               Dsum = zeros(size(parsedData{m,1}{n,1}{1,1}(1,:)));
               increment = 1;
               increment2 = 0;
               ran = 0;
               
               NBins = 0;
               for p=1:length(parsedData{m,1}{n,1}) %Loop over variable names
                   if strcmp(ToLoad{m,3}{p},'Data') % Sum over only data
                       %SummedData = zeros(1,560);
                       SummedData = zeros(1,Options.BinTotal);

                       for q=1:length(DataStructure{m,1}{n,1}.TimeStamp) % Loop over data time
                           if increment+increment2 > length(TimeGrid) % Check if increment is beyond time grid
                               break
                           elseif q==length(DataStructure{m,1}{n,1}.TimeStamp) % Check if increment is last in time grid
                               % Set sumed data bin to running sum
                               SummedData(increment+increment2,:) = Dsum;
                               % Set number of summed bins
                               NBinsInc(increment+increment2,1) = NBins;
                               NewTimeGrid(1,increment+increment2) = TimeGrid(increment+increment2);                        
                           elseif DataStructure{m,1}{n,1}.TimeStamp(q)>TimeGrid(increment+increment2) && increment ==1 && ran ==0 % Check if there is no data before Options.TimeGrid
                               disp('ran DataStructure{m,1}{n,1}.TimeStamp(q)>Options.TimeGrid(increment+increment2) && increment ==1 && ran ==0 % Check if there is no data before Options.TimeGrid')
                               while DataStructure{m,1}{n,1}.TimeStamp(q)>TimeGrid(increment+increment2) % Increment over all missing bins
                                   %SummedData(increment+increment2,:) = NaN(1,560);
                                   SummedData(increment+increment2,:) = NaN(1,Options.BinTotal);
                                   NBinsInc(increment+increment2,1) = NaN;
                                   NewTimeGrid(1,increment+increment2) = TimeGrid(1,increment+increment2);
                                   increment2 = increment2 + 1;
                               end
                               increment2 = increment2-1;
                               increment = increment+1;
                           elseif DataStructure{m,1}{n,1}.TimeStamp(q)<=TimeGrid(increment+increment2) && DataStructure{m,1}{n,1}.TimeStamp(q)>TimeGrid(increment+increment2)-(TimeGrid(2)-TimeGrid(1))% Check if data time increment is inside current time grid section
                               % Add current data bin to sum
                               Dsum = Dsum + parsedData{m,1}{n,1}{p,1}(q,:);
                               % Add number of summed bins
                               NBins = NBins + 1;
                               ran = 1;
                           elseif DataStructure{m,1}{n,1}.TimeStamp(q)>TimeGrid(increment+increment2) % If data time incremten is beyond current time grid section
                               % Set new data bin to sum
                               SummedData(increment+increment2,:) = Dsum;
                               % Set number of summed bins
                               NBinsInc(increment+increment2,1) = NBins;
                               NewTimeGrid(1,increment+increment2) = TimeGrid(increment+increment2);
                               % Reset sum to current bin and NBins
                               Dsum = parsedData{m,1}{n,1}{p,1}(q,:);
                               NBins = 0;
                               % Move to new time grid section
                               increment = increment +1 ;
                           else
                               %do nothing for now
                           end
                       end
                       
                       
                       %Integrate counts in range
                       increment = 1;
                       %clear rangeSummedData
                       %for jj = Options.intRange:Options.intRange:560
                       %for jj = Options.intRange:Options.intRange:Options.BinTotal
                       for jj = 1:1:Options.BinTotal
                           %rangeSummedData(:,increment) = sum(SummedData(:,jj-Options.intRange+1:jj),2);
                           rangeSummedData(:,increment) = sum(SummedData(:,jj-1+1:jj),2);
                           increment = increment+1;
                       end
                       % Take mean of summed data bins and add to structure
                       % cell

                       %DataStructure{m,1}{n,1}.(ToLoad{m,3}{p,:}) = rangeSummedData./NBinsInc;
                       DataStructure{m,1}{n,1}.(ToLoad{m,3}{p,:}) = rangeSummedData;
                       DataStructure{m,1}{n,1}.NBins = NBinsInc;
                       DataStructure{m,1}{n,1}.NewTimeGrid = NewTimeGrid';

                       ToLoad{m,3}{7,1}='NewTimeGrid';%add new varible to ToLoad
                   else
                       % Do nothing for now with other MCS and Power data
                   end
               end        
            else
                % Interpolating 1 dimensional data
                DataStructure{m,1}{n,1} = RecursivelyInterpolateStructure(DataStructure{m,1}{n,1}, ...
                                                                            DataStructure{m,1}{n,1}.TimeStamp,...
                                                                            TimeGrid,...
                                                                            Options.InterpMethod,...
                                                                            Options.Extrapolation);
            end
        end
    end

    %Concatenate days together
    for m=1:size(parsedData,1)
        for n=1:size(parsedData{m,1},1) %Loop over demux options
            for p=1:length(ToLoad{m,3}) %Loop over variable names
                if ~isequal(ToLoad{m,3}{p,1},'Type') && ~isequal(ToLoad{m,3}{p,1},'IsLocked')
                    if strcmp(DataNames{m,1},'MCS') || strcmp(DataNames{m,1},'Power')
                        if i ==1
                            fullDataStructure{m,1}{n,1}.(ToLoad{m,3}{p,1})  = DataStructure{m,1}{n,1}.(ToLoad{m,3}{p,1})';
                        else
                            fullDataStructure{m,1}{n,1}.(ToLoad{m,3}{p,1}) = [fullDataStructure{m,1}{n,1}.(ToLoad{m,3}{p,1}) DataStructure{m,1}{n,1}.(ToLoad{m,3}{p,1})'];
                        end
                    else
                        if i ==1
                            fullDataStructure{m,1}{n,1}.(ToLoad{m,3}{p,1})  = DataStructure{m,1}{n,1}.(ToLoad{m,3}{p,1});
                        else
                            fullDataStructure{m,1}{n,1}.(ToLoad{m,3}{p,1}) = [fullDataStructure{m,1}{n,1}.(ToLoad{m,3}{p,1}) DataStructure{m,1}{n,1}.(ToLoad{m,3}{p,1})];
                        end
                    end
                else
                end
            end
        end
    end

    clear rangeSummedData NBinsInc NewTimeGrid

    end %end big loop over days

%%
% Convert Demux cells to structure
parfor m=1:size(parsedData,1)
    fullDataStructure{m,1} = cell2struct(fullDataStructure{m,1},DemuxNames{m,1});
end
%%
% Pushing all data into a single data structure
DataStructure2 = cell2struct(fullDataStructure,DataNames);

end