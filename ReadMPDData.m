% Written By: Robert Stillwell
% Written For: NCAR
% Modificication Info: Created March 19, 2020
% This function is used to load micropulse dial lidar data.

function [Data] = ReadMPDData(ToLoad,CodePath,DataPath)
%
% Input: ToLoad:   A cell array of strings containing the types of files to
%                  load
%        CodePath: Full path of the location of the processing code
%        DataPath: Full path of the location of the data
%                  
% Output: Data:    Cell array containing all of the requested data
%
%% Loading data
FileCol    = 1;
FileVarCol = 2;
CodeVarCol = 3;
VarTypeCol = 4;
%% 
cd(DataPath)
for m=1:1:size(ToLoad,1) % Looping over filetypes
    s = dir(ToLoad{m,FileCol});   % Finding any availible files
    % Loading all files
    if isempty(s)
        % No files to load
        TimeBounds = linspace(0,23.9,100)';
        % Looping over varaibles
        for p=1:1:size(ToLoad{m,FileVarCol})
            fprintf(['Loading: ',ToLoad{m,FileCol},', default\n'])
            % Pre-allocating data array space
            if strcmp(ToLoad{m,CodeVarCol}{p},'TimeStamp')
                Data{m,1}{p,1} = TimeBounds;
            else
                Data{m,1}{p,1} = TimeBounds.*nan; %#ok<*AGROW>
            end
        end
    else
        % Looping over varaibles
        for p=1:1:size(ToLoad{m,FileVarCol})
            fprintf(['Loading: ',ToLoad{m,FileCol},', ',ToLoad{m,FileVarCol}{p},' as ',ToLoad{m,CodeVarCol}{p},'\n'])
            % Pre-allocating data array space
            Data{m,1}{p,1} = []; 
            for n=1:1:size(s,1)% Looping over availible files
                Filename     = s(n,1).name;
                FileVar      = ToLoad{m,FileVarCol}{p};
                VariableType = ToLoad{m,VarTypeCol}{p};
                % Reading data file
                A = ReadVariable(Filename,FileVar,VariableType);
                % Appending data to the other availible data
                %Data{m,1}{p,1} = [Data{m}{p,1};A];

               % Data{m,1}{p,1} = [Data{m}{p,1}; A];

                sx = size(Data{m}{p,1});
                sy = size(A);
                a = max(sx(2),sy(2));
                Data{m,1}{p,1} = [[Data{m}{p,1} zeros(abs([0 a]-sx))];[A zeros(abs([0 a]-sy))]];
            end
        end
    end
end
cd(CodePath)
end

% This function tries to read a variable of interest. If the type is Double
% it will typecast the variable. Otherwise, if teh type is Native the type
% is maintained from the netcdf file. Finally, if the variable type is
% written as a string, h5read is required because ncread does not recognize
% the string type. 
function [A] = ReadVariable(Filename,FileVar,VariableType)
%
% Inputs: Filename:     String containing the desired file name
%         FileVar:      NetCDF variable name to load
%         VariableType: Data type of variable to read
%
% Outputs: A:           Loaded data
%
%% Checking data inputs 
if nargin == 2
    VariableType = 'Double';
end
%% Loading data
try
    if strcmp(VariableType,'Double')
        A = double(ncread(Filename,FileVar));
    elseif strcmp(VariableType,'Native')
        A = ncread(Filename,FileVar);

%             elseif strcmp(VariableType,'Native')
%         A = double(ncread(Filename,FileVar));
    elseif strcmp(VariableType,'String')
        A = h5read(Filename,['/',FileVar]);
    else
        A = [];
    end
catch
    A = nan;
end

end
