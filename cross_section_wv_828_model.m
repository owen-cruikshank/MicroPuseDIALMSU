function [N_wv,cross_section,lineshape,Line] = cross_section_wv_828_model(T,P,nu_Range,absorption)
%File: cross_section_wv_828_model.m
%Date: 02/14/2021
%Author: Owen Cruikshank
%Inputs:
%   -T:[K] scalar or (range x time) vector of atmopheric temperature as a function
%   of range
%   -P:[atm] scalar or (range x time) vector of atmopheric pressure as a function
%   of range
%   -nu_Range:[1/cm] scaler wavenumber of online laser. Dimentions: (1 x 1 x wavenumber)
%
%Outputs:
%   -sigma: [m^2] The voight absorption cross section of O2 at the input
%   wavenumber arround the 780nm line, dimensions (range x time)
%   -f: [m] Absorption lineshape function, dimensions (range x time)

c = 2.99792458E8;                       %[m/s] speed of light 
kB = 1.38065E-23;                       %[J/K] Boltzman's constant 
h = 6.62607004E-34;                     %[Js] Planck's constant

mwv = 2.991577548987048e-26; % Mass wv molecule [kg]

A = 6.02214e23;
mwv=18.01528/A/1000;


%reference T and P
T0 = 296;                               %[K]
P0 = 1.0;                               %[atm]

parameters = fopen(fullfile('CalibrationData','WV_line_parameters.out'),'r');   %open file containing HITRAN information
fmt = '%1d %1d %f %e %e %f %f %f %f %f';          %format of HITRAN file
WV_parameters =fscanf(parameters,fmt,[10 inf]);     %place HITRAN parameters in vector a
fclose(parameters);                                 %close file
WV_parameters = WV_parameters';                     %transpose matrix to correct format

[rL, tL] = size(T);                                 %length of range vector x length of time vector

lineshape = zeros(rL,tL);
Voight_profile = zeros(rL,tL);
cross_section = zeros(rL,tL);
Line = cell(length(WV_parameters));

nu_Range = nu_Range * 100;                          %change nu_Range from [1/cm] to [1/m]

strength_threshold = 0;                    %[cm / molecule] line strength threshold

t = -10:.2:10;                                      %Relative freqency to integrate over
t = permute(t,[3 1 2]);                             %[none] shift t to put it in third dimestion

increment = 1;

for i = 1:length(WV_parameters)                     %loop over all line parameters

    nu_O2 = WV_parameters(i,3);                     %[1/cm]
    nu_O2 = nu_O2 * 100;                            %[1/m] absoption wavenumber
    S0_O2 = WV_parameters(i,4);                     %[cm/molecule] line strength
    S0_O2 = S0_O2 / 100;                            %[m/molecule] absoption line strength at T=296K
    %S0_O2 = S0_O2 * 0.75;
    gamma_L = WV_parameters(i,6);                   %[1/cm] linewidth
    gamma_L = gamma_L * 100;                        %[1/m] Air (Lorentz) broadened linewidth
    n_air = WV_parameters(i,9);                     %[unitless] linewidth temperature dependence
    E_lower = WV_parameters(i,8);                   %[1/cm] ground state energy
    E_lower = E_lower * 100;                        %[1/m] ground state energy
%     gamma_D = WV_parameters(i,7);                   %[1/cm]self (Doppler) broadened linewidth
%     gamma_D = gamma_D * 100;                        %[1/m]
    delta_air = WV_parameters(i,10);                %[1/cm/atm] pressure shift induce by air, at p=1atm
    delta_air = delta_air * 100;                    %[1/m/atm]
    
    if S0_O2 * 100 > strength_threshold             %Do not compute cross section for line if it id below the threshold

        nuShifted = nu_O2 + delta_air .* P;         %[1/m] (r x t) Shift line center based on atmopheric pressure
        %nuShifted = nu_O2 ;%+ delta_air .* P;         %[1/m] (r x t) Shift line center based on atmopheric pressure
        
        %temperature shifted line strength
        a_o2 = S0_O2.*(T0./T);
        a_o2 = S0_O2.*(T0./T).^1.5;%for WV
        c_o2 = exp(h.*c./kB.*((1./T0)-(1./T)).*E_lower);
        ST_O2 = a_o2.*c_o2;                                 %[m/molecule](t x r) O2 line strength adjusted for temperature shift from T0
        gamma_L_T = gamma_L * (P/P0).*((T0./T).^n_air);     %[1/m](t x r) Lorentz linewidth adjusted for temperature and pressure shift
        gamma_D_T = (nuShifted/c).*sqrt(2*kB*T*log(2)/mwv); %[1/m](t x r) Dopper linewidth due to temperature


        %voight lineshape
        x = ((nu_Range-nuShifted)./gamma_D_T) * sqrt(log(2));   %[none](t x r)
        y = (gamma_L_T./gamma_D_T) * sqrt(log(2));              %[none](t x r)
        K = (ST_O2./gamma_D_T) * sqrt(log(2)/pi);               %[m^2 / molecule](t x r)

        integration_function = exp(-t.^2)./(y.^2 + (x-t).^2);           %[none] create integration function to integrate over

        integralV = trapz(permute(t,[3 2 1]),integration_function,3);   %[none] integrate over t

        f = log(2).*pi^(-3/2).*gamma_L_T./gamma_D_T.^2.*integralV;      %[m] Absorption lineshape 

        Voight = (y/pi).*integralV;                                     %[none](t x r) Voight lineshape

        sigma =  Voight .* K ;                                          %[m^2/molecule](t x r) Voight absorption cross section


        lineshape = f + lineshape;                  %add on to previous lineshape
        Voight_profile = Voight + Voight_profile;   %add on to previous profile
        cross_section = sigma + cross_section;      %add on to previous cross_section
        
        %Save each line parameters
        Line{increment}.lineshape = f;
        Line{increment}.cross_section = sigma;
        Line{increment}.S0=S0_O2;
        Line{increment}.E_lower = E_lower;

        increment = increment+1; %increment line number
    end
end

N_wv = absorption./cross_section; %[molecule/m3] wv number density

end