


% sondeDir = 'C:\Users\Owen\OneDrive - Montana State University\Research\O2 DIAL\Data\NCAR Boulder Data\Soundings\data1OmMj0';
% 
% date_start = datetime(2022,7,27,'TimeZone','UTC');%yyyy,mm,dd
% date_end = datetime(2023,07,21,'TimeZone','UTC');%yyyy,mm,dd
% span_days = date_start:date_end;
% 
% 
% [sonde_datetime,sondeStruc]=m2hatsradiosondeF(sondeDir,span_days)
% 

function [sonde_datetime,sondeStruc] =  m2hatsradiosonde(path,span_days)
%File: COBradiosonde.m
%Date: 07/30/2020
%Author: Owen Cruikshank
%Inputs:
%   -span_days: datetime vector in UTC time of days. ex.: datetime(2020,2,22,'TimeZone','UTC');%yyyy,mm,dd
%   -path: string of file path to radiosonde data folder. Box\Radiosondes\Data\All Data\
%
%Outputs:
%   -sonde_datetime: [datetime] time of radiosonde launch
%   -sondeStruc: Structure of all data recorded by sonde

%File Header:
%    Record name:    Unit:           Data type:          Divisor: Offset:
%    ---------------------------------------------------------------------
%     time            sec             float (4)          1        0       
%     Pscl            ln              short (2)          1        0       
%     T               K               short (2)          10       0       
%     RH              %               short (2)          1        0       
%     v               m/s             short (2)          -100     0       
%     u               m/s             short (2)          -100     0       
%     Height          m               short (2)          1        30000   
%     P               hPa             short (2)          10       0       
%     TD              K               short (2)          10       0       
%     MR              g/kg            short (2)          100      0       
%     DD              dgr             short (2)          1        0       
%     FF              m/s             short (2)          10       0       
%     AZ              dgr             short (2)          1        0       
%     Range           m               short (2)          0.01     0       
%     Lon             dgr             short (2)          100      0       
%     Lat             dgr             short (2)          100      0       
%     SpuKey          bitfield        unsigned short (2) 1        0       
%     UsrKey          bitfield        unsigned short (2) 1        0       
%     RadarH          m               short (2)          1        30000   
% 
% *************************************************************************************************


    sonde_ind = 0;
    
    for day_i = 1:length(span_days)
        
        %oldPath = pwd;
        
        %Finding year folder
        file_year = num2str(year(span_days(day_i)));
        %fullfile(path,file_year,'*.tsv')
        %%%%%fileDir = dir(fullfile(path,file_year,'*.tsv'));
        fileDir = dir(fullfile(path,'*asc.nc'));
        %cd(fullfile(path,file_year))
        %fileDir = dir('*.tsv')
        %cd(oldPath)
        %loop over files in year folder
        
        sonde = cell(length(fileDir),22);
        
        for file_i=1:length(fileDir)
            %COB_2020-07-30-17-2600
            %Find month and day of sonde files
            file_month = str2double(fileDir(file_i).name(30:31));
            file_day = str2double(fileDir(file_i).name(32:33));
            file_year1 = str2double(fileDir(file_i).name(26:29));
            %Check if sondes are within processing day
            if file_month == month(span_days(day_i)) && file_day == day(span_days(day_i)) && file_year1 == year(span_days(day_i))
                %Find hour minute and second of sonde file
                file_hr = str2double(fileDir(file_i).name(35:36));
                file_min = str2double(fileDir(file_i).name(37:38));
                file_sec = str2double(fileDir(file_i).name(39:40));
                %Increment number of sondes found
                sonde_ind = sonde_ind + 1;
                %Create datetime of sonde launch
                sonde_datetime(sonde_ind,1) = datetime(str2double(file_year),file_month,file_day,file_hr,file_min,file_sec,'Timezone','UTC');
                %Read sonde file
                
                sondeFile = fileDir(file_i).name;

                sondeStruc(sonde_ind).time = ncread(fullfile(path,sondeFile),'time'); %time since launch in sec
                sondeStruc(sonde_ind).P = ncread(fullfile(path,sondeFile),'pres');
                sondeStruc(sonde_ind).T = ncread(fullfile(path,sondeFile),'tdry')+273.15;
                sondeStruc(sonde_ind).RH = ncread(fullfile(path,sondeFile),'rh');
                sondeStruc(sonde_ind).Height = ncread(fullfile(path,sondeFile),'alt');
                reference_alt = ncread(fullfile(path,sondeFile),'reference_alt');
                sondeStruc(sonde_ind).Height = sondeStruc(sonde_ind).Height-reference_alt;
                sondeStruc(sonde_ind).vt = ncread(fullfile(path,sondeFile),'vt');
                sondeStruc(sonde_ind).TD = ncread(fullfile(path,sondeFile),'dp')+273.15;

                
                sondeStruc(sonde_ind).e = 0.6108*exp(17.27*(sondeStruc(sonde_ind).TD-273.15)./sondeStruc(sonde_ind).TD);%[kPa]water vapor pressure
                sondeStruc(sonde_ind).AH = 2165*sondeStruc(sonde_ind).e./sondeStruc(sonde_ind).T;%[g/m^3] absolute Humidity 
                sondeStruc(sonde_ind).WV = 6.022e23*sondeStruc(sonde_ind).AH/18.01528;%[1/m^3] WV number density
               
                
%                 sondeStruc(sonde_ind).T
%                 [sondeStruc(sonde_ind).T,outliers(sonde_ind)] = rmoutliers(sondeStruc(sonde_ind).T);
%                 sondeStruc(sonde_ind).P = sondeStruc(sonde_ind).P(outliers(sonde_ind));
%                 sondeStruc(sonde_ind).Height = sondeStruc(sonde_ind).Height(outliers(sonde_ind));
            end
        end
    end
    %If no sonde files are found, set outputs to NaN
    if sonde_ind==0
        sonde_datetime=NaN;
        sondeStruc(1,1:19)=NaN;
    end
end
