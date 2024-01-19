function [HSRL] = HSRL_retrieval_20220909(Counts,Atmosphere,Options,Spectrum)

%this function performs the HSRL retrieval

%load in the calibration constants
%home=pwd;
%cd('Calibration')

    load('CalibrationData\CalibrationScan20220809.mat'); % contains 1 structure called Tables
    % this structure contains the following variables
    %Tables.eta_m (CONSTANT) - molecular efficiency (not used)
    %Tables.eta_c (CONSTANT) - combined efficiency (not used)
    %Tables.Cac (3D TABLE) 
    %Tables.Cmm (3D TABLE)
    %Tables.Cmc (3D TABLE)
    %Tables.Cam (3D TABLE)
    %Tables.Temperature (TABLE DIMENSION 1) [K]
    %Tables.Pressure (TABLE DIMENSION 2) [KPa]
    %Tables.Wavelength (TABLE DIMENSION 3) [m]
  
	calibration_uncertainty = load('CalibrationData\uncertainty_20220809.mat'); % contains 1 structure called calibration_uncertainty
    % this structure contains the following variables
    %calibration_uncertainty.Cac_mean (CONSTANT)
    %calibration_uncertainty.Cac_std (CONSTANT)
    %calibration_uncertainty.Cam_mean (CONSTANT)
    %calibration_uncertainty.Cam_std (CONSTANT)
    %calibration_uncertainty.Cmc_mean (CONSTANT)
    %calibration_uncertainty.Cmc_std (CONSTANT)
    %calibration_uncertainty.Cmm_mean (CONSTANT)
    %calibration_uncertainty.Cmm_std (CONSTANT)
%cd(home);

load('CalibrationData\TransmissionData20220809.mat','Data_Wavelength')

  HSRL.onlineCombinedTransmission = interp1(Data_Wavelength.lambda_on.*10^9,Data_Wavelength.Nc_on,Spectrum.lambda_scan_3D_short);
HSRL.onlineMolecularTransmission = interp1(Data_Wavelength.lambda_on.*10^9,Data_Wavelength.Nm_on,Spectrum.lambda_scan_3D_short);
HSRL.offlineCombinedTransmission = interp1(Data_Wavelength.lambda_off.*10^9,Data_Wavelength.Nc_off,Spectrum.lambda_scan_3D_short_off);
HSRL.offlineMolecularTransmission = interp1(Data_Wavelength.lambda_off.*10^9,Data_Wavelength.Nm_off,Spectrum.lambda_scan_3D_short_off);


for k = 1:length(Options.t_step) % this loop is here because I am running a few different integration times and distances. If you want to use only one, then delete this for loop
    
    %store the pressure and temperature
    [P,T,wavelength] = meshgrid(Tables.Pressure,Tables.Temperature,Tables.Wavelength);

    %store the calibration constants
    Cmm_tab = Tables.Cmm;
    Cmc_tab = Tables.Cmc;
    Cam_tab = Tables.Cam;
    Cac = Tables.Cac(1,1);
    
    % temperature, pressure, and wavelength uncertainty estimate
    sigma_T = 4.6851; %K, uncertainty of the temperature
    dT = 0.001;       %K, small change of temperature for numerical derivative
    sigma_P = 5.1135; % KPa, uncertainty of the pressure
    dP = 0.001;       % KPa, small change of pressure for numerical derivative
    sigma_lambda = 2.3366e-5*1e-9; %m, uncertainty of the wavelength
    dlambda = .000001e-9;% m, small change of wavelength for numerical derivative
    cov_TP = -163.6544; % K-KPa, covariance of the temperature and the pressure
    
    %% find the calibration constants
    Patm = Atmosphere(k).Pressure; % (KPa) define the atmospheric pressure as a variable 
    Tatm = Atmosphere(k).Temperature; % (K) define the atmospheric temperautre as a variable
    lambda = 770.1085e-9*ones(size(Atmosphere(k).Pressure)); % (m) create a 2d variable with a value equal to the target transmitted wavelength
    
    Cmm             = interp3(P,T,wavelength,Cmm_tab,Patm   ,Tatm   ,lambda        ); 
    Cmm_dT_pos      = interp3(P,T,wavelength,Cmm_tab,Patm   ,Tatm+dT,lambda        ); % for numerical derivative wrt temperature
    Cmm_dT_neg      = interp3(P,T,wavelength,Cmm_tab,Patm   ,Tatm-dT,lambda        ); % for numerical derivative wrt temperature
    Cmm_dP_pos      = interp3(P,T,wavelength,Cmm_tab,Patm+dP,Tatm   ,lambda        ); % for numerical derivative wrt pressure
    Cmm_dP_neg      = interp3(P,T,wavelength,Cmm_tab,Patm-dP,Tatm   ,lambda        ); % for numerical derivative wrt pressure
    Cmm_dlambda_pos = interp3(P,T,wavelength,Cmm_tab,Patm   ,Tatm   ,lambda+dlambda); % for numerical derivative wrt wavelength
    Cmm_dlambda_neg = interp3(P,T,wavelength,Cmm_tab,Patm   ,Tatm   ,lambda-dlambda); % for numerical derivative wrt wavelength
    
    Cmc             = interp3(P,T,wavelength,Cmc_tab,Patm,Tatm,lambda);
    Cmc_dT_pos      = interp3(P,T,wavelength,Cmc_tab,Patm   ,Tatm+dT,lambda        ); % for numerical derivative wrt temperature
    Cmc_dT_neg      = interp3(P,T,wavelength,Cmc_tab,Patm   ,Tatm-dT,lambda        ); % for numerical derivative wrt temperature
    Cmc_dP_pos      = interp3(P,T,wavelength,Cmc_tab,Patm+dP,Tatm   ,lambda        ); % for numerical derivative wrt pressure
    Cmc_dP_neg      = interp3(P,T,wavelength,Cmc_tab,Patm-dP,Tatm   ,lambda        ); % for numerical derivative wrt pressure
    Cmc_dlambda_pos = interp3(P,T,wavelength,Cmc_tab,Patm   ,Tatm   ,lambda+dlambda); % for numerical derivative wrt wavelength
    Cmc_dlambda_neg = interp3(P,T,wavelength,Cmc_tab,Patm   ,Tatm   ,lambda-dlambda); % for numerical derivative wrt wavelength
   
    Cam             = interp3(P,T,wavelength,Cam_tab,Patm,Tatm,lambda);
    Cam_dT_pos      = interp3(P,T,wavelength,Cam_tab,Patm   ,Tatm+dT,lambda        ); % for numerical derivative wrt temperature
    Cam_dT_neg      = interp3(P,T,wavelength,Cam_tab,Patm   ,Tatm-dT,lambda        ); % for numerical derivative wrt temperature
    Cam_dP_pos      = interp3(P,T,wavelength,Cam_tab,Patm+dP,Tatm   ,lambda        ); % for numerical derivative wrt pressure
    Cam_dP_neg      = interp3(P,T,wavelength,Cam_tab,Patm-dP,Tatm   ,lambda        ); % for numerical derivative wrt pressure
    Cam_dlambda_pos = interp3(P,T,wavelength,Cam_tab,Patm   ,Tatm   ,lambda+dlambda); % for numerical derivative wrt wavelength
    Cam_dlambda_neg = interp3(P,T,wavelength,Cam_tab,Patm   ,Tatm   ,lambda-dlambda); % for numerical derivative wrt wavelength
    
    clear  Cam_tab Cmm_tab Cmc_tab wavelength % clear the values from the calibration table
    
    % define counts as variables
    Nm_off = Counts(k).Nm_off;
    Nc_off = Counts(k).Nc_off;
    Nm_on = Counts(k).Nm_on;
    Nc_on = Counts(k).Nc_on;

    %% Model the molecular backscatter coefficient and its uncertainty
    %Create Molecular Backscatter Coefficient 
    Beta_mol        = 9.94266e-7*(Patm)   ./(Tatm)   ;
    Beta_mol_dT_pos = 9.94266e-7*(Patm)   ./(Tatm+dT); % for numerical derivative wrt temperature
    Beta_mol_dT_neg = 9.94266e-7*(Patm)   ./(Tatm-dT); % for numerical derivative wrt temperature
    Beta_mol_dP_pos = 9.94266e-7*(Patm+dP)./(Tatm)   ; % for numerical derivative wrt pressure
    Beta_mol_dP_neg = 9.94266e-7*(Patm-dP)./(Tatm)   ; % for numerical derivative wrt pressure
    
    % find the uncertainty for the molecular backscatter coefficient
    d_Beta_mol_d_T = -2*9.94266e-7*(Patm)./(Tatm).^2;
    d_Beta_mol_d_P =    9.94266e-7./(Tatm);
    
    var_Beta_mol = d_Beta_mol_d_T.^2*sigma_T.^2 +d_Beta_mol_d_P.^2*sigma_P.^2 + d_Beta_mol_d_T.*d_Beta_mol_d_P.*cov_TP;
    
    sigma_Beta_mol = abs(sqrt(var_Beta_mol));
    
    %% Backscatter Ratio and its uncertainty
    BR             = 1- (Nm_on .* Cmm             .*Nc_off - Nc_on .* Cmc             .* Nm_off) ./ (Nm_on .* Cam             .*Nc_off - Nc_on .* Nm_off);
    
    BR_dT_pos      = 1- (Nm_on .* Cmm_dT_pos      .*Nc_off - Nc_on .* Cmc_dT_pos      .* Nm_off) ./ (Nm_on .* Cam_dT_pos      .*Nc_off - Nc_on .* Nm_off); % for numerical derivative wrt temperature
    BR_dT_neg      = 1- (Nm_on .* Cmm_dT_neg      .*Nc_off - Nc_on .* Cmc_dT_neg      .* Nm_off) ./ (Nm_on .* Cam_dT_neg      .*Nc_off - Nc_on .* Nm_off); % for numerical derivative wrt temperature
    BR_dP_pos      = 1- (Nm_on .* Cmm_dP_pos      .*Nc_off - Nc_on .* Cmc_dP_pos      .* Nm_off) ./ (Nm_on .* Cam_dP_pos      .*Nc_off - Nc_on .* Nm_off); % for numerical derivative wrt pressure
    BR_dP_neg      = 1- (Nm_on .* Cmm_dP_neg      .*Nc_off - Nc_on .* Cmc_dP_neg      .* Nm_off) ./ (Nm_on .* Cam_dP_neg      .*Nc_off - Nc_on .* Nm_off);  % for numerical derivative wrt pressure
    BR_dlambda_pos = 1- (Nm_on .* Cmm_dlambda_pos .*Nc_off - Nc_on .* Cmc_dlambda_pos .* Nm_off) ./ (Nm_on .* Cam_dlambda_pos .*Nc_off - Nc_on .* Nm_off); % for numerical derivative wrt wavelength
    BR_dlambda_neg = 1- (Nm_on .* Cmm_dlambda_neg .*Nc_off - Nc_on .* Cmc_dlambda_neg .* Nm_off) ./ (Nm_on .* Cam_dlambda_neg .*Nc_off - Nc_on .* Nm_off); % for numerical derivative wrt wavelength
    
    %partial derivatives
    d_BR_d_Nm_off =  (Nc_on  .* Cmc) ./ (Nm_on .* Cam .* Nc_off - Nc_on .* Nm_off) - (1     .* Nc_on  .* (Nm_on .* Cmm .* Nc_off - Nc_on .* Cmc .* Nm_off)) ./ (Nm_on .* Cam .* Nc_off - Nc_on .* Nm_off).^2;
    d_BR_d_Nm_on  = -(Cmm .* Nc_off) ./ (Nm_on .* Cam .* Nc_off - Nc_on .* Nm_off) + (Cam   .* Nc_off .* (Nm_on .* Cmm .* Nc_off - Nc_on .* Cmc .* Nm_off)) ./ (Nm_on .* Cam .* Nc_off - Nc_on .* Nm_off).^2;
    d_BR_d_Nc_off = -(Nm_on  .* Cmm) ./ (Nm_on .* Cam .* Nc_off - Nc_on .* Nm_off) + (Nm_on .* Cam    .* (Nm_on .* Cmm .* Nc_off - Nc_on .* Cmc .* Nm_off)) ./ (Nm_on .* Cam .* Nc_off - Nc_on .* Nm_off).^2;
    d_BR_d_Nc_on  =  (Cmc .* Nm_off) ./ (Nm_on .* Cam .* Nc_off - Nc_on .* Nm_off) - (1     .* Nm_off .* (Nm_on .* Cmm .* Nc_off - Nc_on .* Cmc .* Nm_off)) ./ (Nm_on .* Cam .* Nc_off - Nc_on .* Nm_off).^2;
    d_BR_d_Cmm    = -(Nm_on .* Nc_off) ./ (Nm_on .* Cam .* Nc_off - Nc_on .* Nm_off);
    d_BR_d_Cmc    =  (Nc_on .* Nm_off) ./ (Nm_on .* Cam .* Nc_off - Nc_on .* Nm_off);
    d_BR_d_Cam    = (Nm_on .* Nc_off .* (Nm_on .* Cmm .* Nc_off - Nc_on .* Cmc .* Nm_off)) ./ (Nm_on .* Cam .* Nc_off - Nc_on .* Nm_off).^2;
    d_BR_d_T      = (BR_dT_pos      - BR_dT_neg     ) ./ (2*dT);
    d_BR_d_P      = (BR_dP_pos      - BR_dP_neg     ) ./ (2*dP);
    d_BR_d_lambda = (BR_dlambda_pos - BR_dlambda_neg) ./ (2*dlambda);

    % uncertainties
    sigma_Nm_off = Counts(k).sigma_Nm_off;
    sigma_Nm_on  = Counts(k).sigma_Nm_on;
    sigma_Nc_off = Counts(k).sigma_Nc_off;
    sigma_Nc_on  = Counts(k).sigma_Nc_on;
    sigma_Cmm    = calibration_uncertainty.Cmm_std;
    sigma_Cmc    = calibration_uncertainty.Cmc_std;
    sigma_Cam    = calibration_uncertainty.Cam_std;
    
    % combining terms to find the variance
    var_BR = (d_BR_d_Nm_off) .^2 .* (sigma_Nm_off) .^2 ...
           + (d_BR_d_Nm_on ) .^2 .* (sigma_Nm_on ) .^2 ...   
           + (d_BR_d_Nc_off) .^2 .* (sigma_Nc_off) .^2 ...    
           + (d_BR_d_Nc_on ) .^2 .* (sigma_Nc_on ) .^2 ...    
           + (d_BR_d_Cmm   ) .^2 .* (sigma_Cmm   ) .^2 ...    
           + (d_BR_d_Cmc   ) .^2 .* (sigma_Cmc   ) .^2 ...    
           + (d_BR_d_Cam   ) .^2 .* (sigma_Cam   ) .^2 ...    
           + (d_BR_d_lambda) .^2 .* (sigma_lambda) .^2 ...    
           + (d_BR_d_T     ) .^2 .* (sigma_T     ) .^2 ...    
           + (d_BR_d_P     ) .^2 .* (sigma_P     ) .^2 ...
           + (d_BR_d_T) .* (d_BR_d_P) .* cov_TP;
    % uncertainty
    sigma_BR = abs(sqrt(var_BR));

%% find the aerosol backscatter coefficient
    Beta_aer = Beta_mol .* (BR-1);
    
    Beta_aer_dT_pos      = Beta_mol_dT_pos .* (BR_dT_pos     -1); % for numerical derivative wrt temperature
    Beta_aer_dT_neg      = Beta_mol_dT_neg .* (BR_dT_neg     -1); % for numerical derivative wrt temperature
    Beta_aer_dP_pos      = Beta_mol_dP_pos .* (BR_dP_pos     -1); % for numerical derivative wrt pressure
    Beta_aer_dP_neg      = Beta_mol_dP_neg .* (BR_dP_neg     -1); % for numerical derivative wrt pressure
    Beta_aer_dlambda_pos = Beta_mol        .* (BR_dlambda_pos-1); % for numerical derivative wrt wavelength
    Beta_aer_dlambda_neg = Beta_mol        .* (BR_dlambda_neg-1); % for numerical derivative wrt wavelength

    %partial derivatives
    d_Beta_aer_d_Nm_off = Beta_mol .* d_BR_d_Nm_off;
    d_Beta_aer_d_Nm_on  = Beta_mol .* d_BR_d_Nm_on;
    d_Beta_aer_d_Nc_off = Beta_mol .* d_BR_d_Nc_off;
    d_Beta_aer_d_Nc_on  = Beta_mol .* d_BR_d_Nc_on;
    d_Beta_aer_d_Cmm    = Beta_mol .* d_BR_d_Cmm;
    d_Beta_aer_d_Cmc    = Beta_mol .* d_BR_d_Cmc;
    d_Beta_aer_d_Cam    = Beta_mol .* d_BR_d_Cam;
    d_Beta_aer_d_T      = (Beta_aer_dT_pos      - Beta_aer_dT_neg     ) ./ (2*dT);
    d_Beta_aer_d_P      = (Beta_aer_dP_pos      - Beta_aer_dP_neg     ) ./ (2*dP);
    d_Beta_aer_d_lambda = (Beta_aer_dlambda_pos - Beta_aer_dlambda_neg) ./ (2*dlambda);
    
    %combine terms to find the variance
    var_Beta_aer = (d_Beta_aer_d_Nm_off) .^2 .* (sigma_Nm_off) .^2 ...
                 + (d_Beta_aer_d_Nm_on ) .^2 .* (sigma_Nm_on ) .^2 ...   
                 + (d_Beta_aer_d_Nc_off) .^2 .* (sigma_Nc_off) .^2 ...    
                 + (d_Beta_aer_d_Nc_on ) .^2 .* (sigma_Nc_on ) .^2 ...    
                 + (d_Beta_aer_d_Cmm   ) .^2 .* (sigma_Cmm   ) .^2 ...    
                 + (d_Beta_aer_d_Cmc   ) .^2 .* (sigma_Cmc   ) .^2 ...    
                 + (d_Beta_aer_d_Cam   ) .^2 .* (sigma_Cam   ) .^2 ...    
                 + (d_Beta_aer_d_lambda) .^2 .* (sigma_lambda) .^2 ...    
                 + (d_Beta_aer_d_T     ) .^2 .* (sigma_T     ) .^2 ...    
                 + (d_Beta_aer_d_P     ) .^2 .* (sigma_P     ) .^2 ...
                 + (d_Beta_aer_d_T) .* (d_Beta_aer_d_P) .* cov_TP;
       
    % determine the uncertainty         
    sigma_Beta_aer = abs(sqrt(var_Beta_aer));
    
	%% find the molecular counts 
    % this is for finding the aerosol optical depth
    Nmol =(1./(Cmc.*Cam-Cmm)).*(Cam.*Nc_off - Nc_on./Nm_on.*Nm_off);
     
    %% store the variables
    %HSRL is a stucture containing my primary HSRL data products
    HSRL(k).Bm = Beta_mol; % molecular backscatter coefficient
    HSRL(k).Ba  = Beta_aer; % aerosol backscatter coefficient
    HSRL(k).BSR        = BR; % backscatter ratio
    HSRL(k).sigma_Beta_mol = sigma_Beta_mol; % uncertainty
    HSRL(k).sigma_Beta_aer = sigma_Beta_aer; % uncertainty
    HSRL(k).sigma_BR       = sigma_BR; % uncertainty
   % HSRL(k).time  = Counts(k).time; % save time as a dimension
   % HSRL(k).range = Counts(k).range; % save range as a dimension
%    HSRL(k).BR_aer = HSRL(k).BR-1; % aerosol-to-molecular backscatter ratio
    HSRL(k).Cac = Cac; % calibration constant
    HSRL(k).Cam = Cam; % calibration constant
    HSRL(k).Cmc = Cmc; % calibration constant
    HSRL(k).Cmm = Cmm; % calibration constant
    HSRL(k).eta_m = Tables.eta_m;
    HSRL(k).eta_c = Tables.eta_c;
    HSRL(k).Nmol = Nmol; % molecular counts
end
end