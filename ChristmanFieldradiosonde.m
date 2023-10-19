
function [sonde_datetime,sondeStruc] =  ChristmanFieldradiosonde(path,span_days)
%File: ChristmanFieldradiosonde.m
%Date: 12/20/21
%Author: Owen Cruikshank
%Inputs:
%   -span_days: datetime vector in UTC time of days. ex.: datetime(2020,2,22,'TimeZone','UTC');%yyyy,mm,dd
%   -path: string of file path to radiosonde data folder. Box\Radiosondes\Data\All Data\
%
%Outputs:
%   -sonde_datetime: [datetime] time of radiosonde launch
%   -sondeStruc: Structure of all data recorded by sonde

%File Header:
% Station name                                 	PRE-CIP_ChristmanField
% System trademark and model                   	MW41
% Sonde type                                   	RS41-SGP
% Sonde serial number                          	S0651213
% Ground check device                          	RI41
% Balloon release date and time                	2021-06-20T17:58:39
% Release point latitude                       	40.590000°N
% Release point longitude                      	105.141500°W
% Release point height from sea level          	1571.9 m
% Surface pressure                             	840.4 hPa
% Surface temperature                          	23.6 °C
% Surface relative humidity                    	32 %
% Surface wind speed                           	2.7 m/s
% Surface wind direction                       	138°
% Surface pressure                             	840.4 hPa
% Pref pressure                                	839.95 hPa
% P correction (Pref - P)                      	1.45 hPa
% Tref temperature                             	/////
% Uref humidity                                	Desiccant 0 %Rh
% U correction (Uref1 - U1)                    	0.5 %Rh
% Average ascent rate                          	4.3 m/s
% Terminating altitude                         	4980 m
% Sounding length                              	00:13:18 hh:mm:ss
% Reason for termination                       	WeakOrFadingSignal
% Reason for sounding failure                  	/////
% Height and pressure in messages is based on  	P Sensor
% Software version                             	MW41 2.17.0
% 
% Elapsed time HeightMSL     P Temp  Dewp   RH Speed   Dir AscRate       Lat         Lon GpsHeightMSL
%            s         m   hPa   °C    °C    %   m/s     °     m/s         °           °            m
    sonde_ind = 0;
    
    for day_i = 1:length(span_days)
        
        %oldPath = pwd;
        
        %Finding year folder
        file_year = num2str(year(span_days(day_i)));
        %fullfile(path,file_year,'*.tsv')
        fileDir = dir(fullfile(path,'*.txt'));
        %cd(fullfile(path,file_year))
        %fileDir = dir('*.tsv')
        %cd(oldPath)
        %loop over files in year folder
        
        sonde = cell(length(fileDir),15);
        
        for file_i=1:length(fileDir)
            %edt_20210620_1758.txt
            %Find month and day of sonde files
            file_month = str2double(fileDir(file_i).name(9:10));
            file_day = str2double(fileDir(file_i).name(11:12));
            %Check if sondes are within processing day
            if file_month == month(span_days(day_i)) && file_day == day(span_days(day_i))
                %Find hour minute and second of sonde file
                file_hr = str2double(fileDir(file_i).name(14:15));
                file_min = str2double(fileDir(file_i).name(16:17));
                %%file_sec = str2double(fileDir(file_i).name(21:22));
                file_sec = 30;
                %Increment number of sondes found
                sonde_ind = sonde_ind + 1;
                %Create datetime of sonde launch
                sonde_datetime(sonde_ind,1) = datetime(str2double(file_year),file_month,file_day,file_hr,file_min,file_sec,'Timezone','UTC');
                %Read sonde file
                formatSpec='%f %f %f %f %f %f %f %f %f %f %f %f';
                sondeFile = fileread(fullfile(fileDir(file_i).folder,fileDir(file_i).name));
                sonde(sonde_ind,1:end-3) = textscan(sondeFile,formatSpec,'Delimiter',' ','MultipleDelimsAsOne',1,'HeaderLines',30);
                
                %Create sonde structure with headings
                %%headings = {'time','Pscl','T','RH','v','u','Height','P','TD','MR','DD','FF','AZ','Range','Lon','Lat','SpuKey','UsrKey','RadarH','e','AH','WV'};
                headings = {'time','Height','P','T','TD','RH','Speed','Dir','AscRate','Lat','Lon','GPSHeightMSL','e','AH','WV'};
                sondeStruc(sonde_ind) = cell2struct(sonde(sonde_ind,:),headings,2);
                
                notOutliers = find(sondeStruc(sonde_ind).T~=-32768);
                
                sondeStruc(sonde_ind).T = sondeStruc(sonde_ind).T(notOutliers)+273.15;
                sondeStruc(sonde_ind).P = sondeStruc(sonde_ind).P(notOutliers);
                sondeStruc(sonde_ind).Height = sondeStruc(sonde_ind).Height(notOutliers);
                sondeStruc(sonde_ind).Height = sondeStruc(sonde_ind).GPSHeightMSL(notOutliers);
                %sondeStruc(sonde_ind).Height = sondeStruc(sonde_ind).Height(notOutliers)+139+400;
                
                sondeStruc(sonde_ind).TD = sondeStruc(sonde_ind).TD(notOutliers)+273.15;
                sondeStruc(sonde_ind).RH = sondeStruc(sonde_ind).RH(notOutliers);
                
                
                sondeStruc(sonde_ind).e = 0.6108*exp(17.27*(sondeStruc(sonde_ind).TD-273.15)./sondeStruc(sonde_ind).TD);%[kPa]water vapor pressure
                sondeStruc(sonde_ind).AH = 2165*sondeStruc(sonde_ind).e./sondeStruc(sonde_ind).T;%[g/m^3] absolute Humidity 
                sondeStruc(sonde_ind).WV = 6.022e23*sondeStruc(sonde_ind).AH/18.01528;%[1/m^3] WV number density
               
                f = 1.0016+3.15.*10.^-6.*sondeStruc(sonde_ind).P-.074.*sondeStruc(sonde_ind).P.^-1;
                ew =f.*6.1094.*exp((17.625.*(sondeStruc(sonde_ind).T-273.15)./(243.04+(sondeStruc(sonde_ind).T-273.15))));
                ew =6.1094.*exp((17.625.*(sondeStruc(sonde_ind).T-273.15)./(243.04+(sondeStruc(sonde_ind).T-273.15))));
                ew = (1.0007+3.46e-6.*sondeStruc(sonde_ind).P).*6.1121.*exp(17.502.*(sondeStruc(sonde_ind).T-273.15)./(240.97+(sondeStruc(sonde_ind).T-273.15)));
                e=ew.*sondeStruc(sonde_ind).RH./100;
                sondeStruc(sonde_ind).AH = e./461.5./(sondeStruc(sonde_ind).T-273.15); %g/m^3


                R = 8.31446261815324;%J/K/mol, m^3Pa/K/mol

                sondeStruc(sonde_ind).AH = (e./10).*18.01528./R./sondeStruc(sonde_ind).T*1000; %g/m^3

                A = 6.02214e23; %avagadro number
                mwv=18.01528/A; %g/molecule
                sondeStruc(sonde_ind).WV = sondeStruc(sonde_ind).AH./mwv; %1/m^3
                
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
