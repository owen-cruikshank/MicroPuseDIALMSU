function [absorption,cross_section,lineshape,Line] = absorption_O2_770_model(T,P,nu_Range,WV)
%File: absorption_O2_770_model.m
%Date: 02/28/2020
%Author: Owen Cruikshank
%Inputs:
%   -T:[K] scalar or (range x time) vector of atmopheric temperature as a function
%   of range
%   -P:[atm] scalar or (range x time) vector of atmopheric pressure as a function
%   of range
%   -nu_Range:[1/cm] scaler wavenumber of online laser. Dimentions: (1 x 1 x wavenumber)
%   -WV:[1/m^3] water vapor number density. Dimentions: (range x time)
%
%Outputs:
%   -absorption: [m-1] the atmopheric absorption of O2 at the input
%   wavenumber arround the 780nm line, dimensions (range x time)
%   -sigma: [m^2] The voight absorption cross section of O2 at the input
%   wavenumber arround the 780nm line, dimensions (range x time)
%   -f: [m] Absorption lineshape function, dimensions (range x time)

c = 2.99792458E8;                       %[m/s] speed of light 
kB = 1.38065E-23;                       %[J/K] Boltzman's constant 
h = 6.62607004E-34;                     %[Js] Planck's constant

q_O2 = .2095;                           %[unitless] O2 atmospheric mixing ratio 

mo2 = 5.313525282632483E-26;            % Mass O2 molecule [kg]

%reference T and P
T0 = 296;                               %[K]
P0 = 1.0;                               %[atm]

% parameters = fopen(fullfile('CalibrationData','O2_line_parameters.out'),'r');   %open file containing HITRAN information
% fmt = '%1d %1d %f %e %e %f %f %f %f %f';          %format of HITRAN file
% O2_parameters =fscanf(parameters,fmt,[10 inf]);     %place HITRAN parameters in vector a
% fclose(parameters);                                 %close file
% O2_parameters = O2_parameters';                     %transpose matrix to correct format

 O2_parameters=[7	1	12990.4577700000	4.88900000000000e-26	0.0219200000000000	0.0312000000000000	0.0340000000000000	1420.76310000000	0.630000000000000	-0.00930000000000000;
 7	1	12990.5018320000	3.74900000000000e-27	0.0171900000000000	0.0491000000000000	0.0490000000000000	1635.06590000000	0.740000000000000	-0.00730000000000000];


%2016
O2_parameters=[
7 1	12990.45777	 4.889e-26 0.02192 0.0312 0.034	1420.7631 0.63 -0.0093;
7 1	12990.501832 3.749e-27 0.01719 0.0491 0.049	1635.0659 0.74 -0.0073;
7 1 12988.722523 5.059e-26 0.02287 0.0312 0.034 1422.4983 0.63 -0.0093;
7 1 12984.267206 2.674e-27 0.01352 0.0507 0.052 1600.1301 0.73 -0.0061;
7 1 12986.261218 3.404e-27 0.01705 0.0507 0.052 1598.1361 0.73 -0.0071];

 %2012
O2_parameters=[
7 1 12990.457779 4.86e-26 0.02192 0.0312 0.034 1420.766 0.63 -0.0093;
7 1 12990.502232 3.749e-27 0.0172 0.0491 0.049 1635.0686 0.74 -0.0073;
7 1 12988.722531 5.029e-26 0.02287 0.0312 0.034 1422.5012 0.63 -0.0093;
7 1 12984.267824 2.674e-27 0.01352 0.0507 0.052 1600.1329 0.73 -0.0061;
7 1 12986.261834 3.404e-27 0.01705 0.0507 0.052 1598.1388 0.73 -0.0071];

f = fopen('5d41b591.par','r');
formatSpec = '%1d%1d%f %e %e%6f%4f %f%4f%8f       b      %1f       X      %1f                %1c %2d%1c %2d     %1c%14c %d %d %f %f';
O2_parameters = fscanf(f,formatSpec,[35 Inf]);
O2_parameters = O2_parameters';
fclose(f);

O2_parameters = O2_parameters([3 4 8 11 12],:);



% O2_parameters=[7	1	12990.4577700000	4.88900000000000e-26	0.0219200000000000	0.0312000000000000	0.0340000000000000	1420.76310000000	0.630000000000000	-0.00930000000000000]

% parameters = fopen(fullfile('CalibrationData','O2_line_parameters.out'),'r');   %open file containing HITRAN information
% fmt = '%1d %1d %f %e %e %f %f %f %f %f';          %format of HITRAN file
% O2_parameters =fscanf(parameters,fmt,[10 inf]);     %place HITRAN parameters in vector a
% fclose(parameters);                                 %close file
% O2_parameters = O2_parameters';                     %transpose matrix to correct format


[rL, tL] = size(T);                                 %length of range vector x length of time vector

lineshape = zeros(rL,tL);
cross_section = zeros(rL,tL);

nu_Range = nu_Range * 100;                          %change nu_Range from [1/cm] to [1/m]

strength_threshold = 1*10^(-26);                    %[cm / molecule] line strength threshold
strength_threshold = 1*10^(-25);
strength_threshold = 3*10^(-26);
strength_threshold = 0;

t = -10:.2:10;                                      %Relative freqency to integrate over
t = permute(t,[3 1 2]);                             %[none] shift t to put it in third dimestion
t1 =permute(t,[3 2 1]);

%Preallocate
% linesThreshold=0;
% for ii = 1:length(O2_parameters) 
%     if  O2_parameters(ii,4)> strength_threshold
%         linesThreshold=linesThreshold+1;
%     end
% end
% Line = cell(1,linesThreshold);
Line = cell(1,2);

increment = 1;
for i = 1:size(O2_parameters,1)                     %loop over all line parameters

    nu_O2 = O2_parameters(i,3);                     %[1/cm]
    nu_O2 = nu_O2 * 100;                            %[1/m] absoption wavenumber
    S0_O2 = O2_parameters(i,4);                     %[cm/molecule] line strength
    S0_O2 = S0_O2 / 100;                            %[m/molecule] absoption line strength at T=296K
    gamma_L = O2_parameters(i,6);                   %[1/cm] linewidth
    gamma_L = gamma_L * 100;                        %[1/m] Air (Lorentz) broadened linewidth
    n_air = O2_parameters(i,9);                     %[unitless] linewidth temperature dependence
    E_lower = O2_parameters(i,8);                   %[1/cm] ground state energy
    E_lower = E_lower * 100;                        %[1/m] ground state energy
    delta_air = O2_parameters(i,10);                %[1/cm/atm] pressure shift induce by air, at p=1atm
    delta_air = delta_air * 100;                    %[1/m/atm]
    
    if S0_O2 * 100 > strength_threshold             %Do not compute cross section for line if it id below the threshold

        
        nuShifted = nu_O2 + delta_air .* P;         %[1/m] (r x t) Shift line center based on atmopheric pressure
      
        %temperature shifted line strength
        %ST_O2 = S0_O2.*(T0./T).*exp(h.*c./kB.*((1./T0)-(1./T)).*E_lower);                                 %[m/molecule](t x r) O2 line strength adjusted for temperature shift from T0
        ST_O2 = S0_O2.*(T0./T).*exp(h.*c./kB.*((1./T0)-(1./T)).*E_lower).*((1-exp(-h*c*nu_O2./kB./T))./(1-exp(-h*c*nu_O2./kB./T0))); 

        % Q296 = TIPS2017(O2_parameters(i,1),O2_parameters(i,2),T0);
        % Q = TIPS2017(O2_parameters(i,1),O2_parameters(i,2),T);
        % ST_O2 = S0_O2.*(Q296./Q).*exp(h.*c./kB.*((1./T0)-(1./T)).*E_lower).*((1-exp(-h*c*nu_O2./kB./T))./(1-exp(-h*c*nu_O2./kB./T0)));
        
        gamma_L_T = gamma_L * (P/P0).*((T0./T).^n_air);     %[1/m](t x r) Lorentz linewidth adjusted for temperature and pressure shift
        gamma_D_T = (nuShifted/c).*sqrt(2*kB*T*log(2)/mo2); %[1/m](t x r) Dopper linewidth due to temperature

        %voight lineshape
        x = ((nu_Range-nuShifted)./gamma_D_T) * sqrt(log(2));   %[none](t x r)
        y = (gamma_L_T./gamma_D_T) * sqrt(log(2));              %[none](t x r)
        K = (ST_O2./gamma_D_T) * sqrt(log(2)/pi);               %[m^2 / molecule](t x r)

        integration_function = exp(-t.^2)./(y.^2 + (x-t).^2);           %[none] create integration function to integrate over

        integralV = trapz(t1,integration_function,3);   %[none] integrate over t

        f = log(2).*pi^(-3/2).*gamma_L_T./gamma_D_T.^2.*integralV;      %[m] Absorption lineshape 

        Voight = (y/pi).*integralV;                                     %[none](t x r) Voight lineshape

        lineshape = f + lineshape;                  %add on to previous lineshape
       cross_section = Voight .* K + cross_section;      %add on to previous cross_section
        
        %Save each line parameters
        %N_o2 = ((P*101325)./(kB*T)-WV) * q_O2; %[molecule/m^3](t x r)O2 number density from atmopsheric number density and O2 mixing ratio
        %Line{increment}.absorption = Voight .* K .*N_o2;
        Line{increment}.lineshape = f;
        Line{increment}.cross_section = Voight .* K;
        Line{increment}.S0=S0_O2;
        Line{increment}.E_lower = E_lower;
        
        Line{increment}.istotope = O2_parameters(i,2); 
        Line{increment}.v = O2_parameters(i,12);
        Line{increment}.quanta = [char(O2_parameters(i,13)) num2str(O2_parameters(i,14)) char(O2_parameters(i,15)) num2str(O2_parameters(i,14))];

        Line{increment}.a =  Voight .* K.* ((P*101325)./(kB*T)-WV) * q_O2; %[1/m](t x r)absorption coefficeint of oxygen in the atmosphere at specificed wavenumber

        increment = increment+1; %increment line number
    end
end

absorption = cross_section .* ((P*101325)./(kB*T)-WV) * q_O2; %[1/m](t x r)absorption coefficeint of oxygen in the atmosphere at specificed wavenumber

end