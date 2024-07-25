function [absorption,cross_section,lineshape,Line] = absorption_O2_770_modelSD(T,P,nu_Range,WV,Constant)
%File: absorption_O2_770_modelSD.m
%Date: 06/01/2024
%Author: Owen Cruikshank
%Description: Calculate Speed Dependent Voight profile
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

c = Constant.c;
kB = Constant.kB;
h = Constant.h;
q_O2 = Constant.q_O2;
mo2 = Constant.mo2;

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

f = fopen(fullfile('CalibrationData','5d41b591.par'),'r');
formatSpec = '%1d%1d%f %e %e%6f%4f %f%4f%8f       b      %1f       X      %1f                %1c %2d%1c %2d     %1c%14c %d %d %f %f';
O2_parameters = fscanf(f,formatSpec,[35 Inf]);
O2_parameters = O2_parameters';
fclose(f);

%O2_parameters = O2_parameters([3 4 8 11 12],:);

O2_parameters = O2_parameters(11,:);



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

    nu_O2 = 12990.45772*100;
    S0_O2 = 4.889e-26/100;
    gamma_L = 0.03223*100;
    n_air = 0.630;
    E_lower = 1420.7631*100;
    delta_air = -0.01066*100;

    
    if S0_O2 * 100 > strength_threshold             %Do not compute cross section for line if it id below the threshold

        
        nuShifted = nu_O2 + delta_air .* P;         %[1/m] (r x t) Shift line center based on atmopheric pressure
      
        %temperature shifted line strength
        %ST_O2 = S0_O2.*(T0./T).*exp(h.*c./kB.*((1./T0)-(1./T)).*E_lower);                                 %[m/molecule](t x r) O2 line strength adjusted for temperature shift from T0
        ST_O2 = S0_O2.*(T0./T).*exp(h.*c./kB.*((1./T0)-(1./T)).*E_lower).*((1-exp(-h*c*nu_O2./kB./T))./(1-exp(-h*c*nu_O2./kB./T0))); 

         % Q296 = TIPS2017(O2_parameters(i,1),O2_parameters(i,2),T0);
         % Q = TIPS2017(O2_parameters(i,1),O2_parameters(i,2),T);
         % ST_O2 = S0_O2.*(Q296./Q).*exp(h.*c./kB.*((1./T0)-(1./T)).*E_lower).*((1-exp(-h*c*nu_O2./kB./T))./(1-exp(-h*c*nu_O2./kB./T0)));

        gamma_L_T = gamma_L * (P/P0).*((T0./T).^n_air);     %[1/m](t x r) Lorentz linewidth adjusted for temperature and pressure shift
        Pwv = 0.01*P;
        Pwv = 0;
        gamma_Lwv = 0.0343;
        n_wv = 0.61;
        gamma_L_Twv = gamma_L * ((P-Pwv)/P0).*((T0./T).^n_air)+gamma_Lwv * ((Pwv)/P0).*((T0./T).^n_wv);
        gamma_D_T = (nuShifted/c).*sqrt(2*kB*T*log(2)/mo2); %[1/m](t x r) Dopper linewidth due to temperature

        %voight lineshape
        x = ((nu_Range-nuShifted)./gamma_D_T) * sqrt(log(2));   %[none](t x r)
        y = (gamma_L_T./gamma_D_T) * sqrt(log(2));              %[none](t x r)
        K = (ST_O2./gamma_D_T) * sqrt(log(2)/pi);               %[m^2 / molecule](t x r)

        integration_function = exp(-t.^2)./(y.^2 + (x-t).^2);           %[none] create integration function to integrate over

        integralV = trapz(t1,integration_function,3);   %[none] integrate over t

        f = log(2).*pi^(-3/2).*gamma_L_T./gamma_D_T.^2.*integralV;      %[m] Absorption lineshape 

        Voight = (y/pi).*integralV;                                     %[none](t x r) Voight lineshape


        %%
        %K = real(erf(x+1i*y));
        %%%K = real(w(x+1i*y));
%        K1 = real(exp(-(x+1i*y).^2).*erf(-1i.*(x+1i*y)));
        %%%gv=sqrt(log(2)/pi)./gamma_D_T.*K;
        %%%f=gv;
        %%
        %SD voight
        S = 0.1;
        gamma2 = S;
        delta = 1/(4*log(2)) .* (gamma_D_T./gamma2).^2;
        alpha = 2*y*sqrt(delta)-3/2;
        beta = 2*x*sqrt(delta);
        zplus = sqrt(delta) + 1/sqrt(2).*sqrt(sqrt((alpha+delta).^2+beta.^2)+alpha+delta) +(1i.*beta./sqrt(2))./(sqrt(sqrt((alpha+delta).^2+beta.^2)+alpha+delta));
        zminus = -sqrt(delta) + 1/sqrt(2).*sqrt(sqrt((alpha+delta).^2+beta.^2)+alpha+delta) +(1i.*beta./sqrt(2))./(sqrt(sqrt((alpha+delta).^2+beta.^2)+alpha+delta));
        %Q = real(erf(1i.*zminus)-erf(1i.*zplus));
        Q = real(w(1i.*zminus)-w(1i.*zplus));
        g = sqrt(log(2)/pi).*Q./gamma_D_T;
        %%
        %HT
        %SDV
        eta = 0; %partial correlation between velocity and rotational state changes due to collisions.
        Delta_2 = 0; %quadratic speed dependence 
        gamma_2 = 0; %quadratic relaxation

        Delta_2 = 0.0001; %quadratic speed dependence 
        gamma_2 = 0.0001; %quadratic relaxation

        Delta_2 = 5e-5; %quadratic speed dependence 
        gamma_2 = 0; %quadratic relaxation
        nu_VC = 0;%Nelkin–Ghatak hard collision model
        gamma_0 = gamma_L_T;%lorentzian half width
        Delta_0 = delta_air.*P;
        nu_0 = nu_O2;
        nu = nu_Range;

        nu_a0 = c./(sqrt(log(2)).*nu_0).*gamma_D_T;%averagegaussian speed
        C_0 = gamma_0 + 1i.*Delta_0;
        C_2 = gamma_2 + 1i.*Delta_2;

        C_0t = (1-eta).*(C_0-3*C_2/2)+nu_VC;
        C_2t = (1-eta).*C_2;
        X = (-1i.*(nu_0-nu)+C_0t)./C_2t;
        Y = (nu_0.*nu_a0./(2*c*C_2)).^2;
        zplus = sqrt(X+Y)+sqrt(Y);
        zminus = sqrt(X+Y)-sqrt(Y);
        A = sqrt(pi)*c./(nu_0.*nu_a0) .* (w(1i.*zminus)-w(1i.*zplus));
        B = nu_a0.^2./C_2t .* (-1 + sqrt(pi)./(2.*sqrt(Y)).*(1-zminus.^2).*w(1i.*zminus) - sqrt(pi)./(2.*sqrt(Y)).*(1-zplus.^2).*w(1i.*zplus));
        
        f_HTP = (1/pi).*real(A ./(1-(nu_VC-eta.*(C_0t-3.*C_2t/2).*A+(eta.*C_2t./nu_a0.^2).*B) ));



        %%
        %SD Rautian Drouin et al
        v = linspace(-4,4,17);
        %v = linspace(-40,40,100);
        xf = 1;
        xs = 0;
        df = -0.01066*100;
        dfp = 5e-5*100;
        ds = -.00660*100;
        dsp = 3e-5*100;
        gammaf = 0.03223*100;
        gammas = 0.03356*100;
        nf = 0.630;
        ns = 0.750;
        H =0;
        S = 0.1;
        S=0;

        vMostPobalby = sqrt(2*kB*T/mo2);

        %v = v./vMostPobalby;

        x = (nu-nu_0-P.*(xf.*(df+(T-T0).*dfp)+xs.*(ds+(T-T0).*dsp))).*sqrt(log(2))./gamma_D_T ;
        y = (P.*(xf.*gammaf.*(T0./T).^nf+xs.*gammas.*(T0./T).^ns)).*sqrt(log(2))./gamma_D_T   ;

        x = (nu-nu_0-P.*(xf.*(df+(T-T0).*dfp)+xs.*(ds+(T-T0).*dsp)))./gamma_D_T ;
        y = (P.*(xf.*gammaf.*(T0./T).^nf+xs.*gammas.*(T0./T).^ns))./gamma_D_T   ;
        F = (2./pi.^(3/2)).*trapz(v,v.*exp(-v.^2).*atan((x+v)./(y.*(1+S.*(v.^2-3/2))+H)));

        c1=-2.658e4;
        c2=3.36e6;
        muso = 0.027430;
        GEVoverQa = 1.1973e28;%cm^2

        m = -31;

        %Gev = 8*pi^3*muso^2./(3*Constant.h*Constant.c*Qevrs)
        %I = Gev.*nu_0.*SHL.*(1+c1*m+c2*m^2).*exp(-E_lower/Constant.kB./T)
        I = ST_O2;
        k = q_O2 .* F .* I.*((P*Constant.ATMtoPA)./(kB*T)-WV);


        %%
        %Voight = g;

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

        Line{increment}.a =  Voight .* K.* ((P*Constant.ATMtoPA)./(kB*T)-WV) * q_O2; %[1/m](t x r)absorption coefficeint of oxygen in the atmosphere at specificed wavenumber

        increment = increment+1; %increment line number
    end
end

absorption = cross_section .* ((P*Constant.ATMtoPA)./(kB*T)-WV) * q_O2; %[1/m](t x r)absorption coefficeint of oxygen in the atmosphere at specificed wavenumber

end

function  erf = w(z)
    t = -10:.2:10;                                      %Relative freqency to integrate over
    t = permute(t,[3 1 2]);                             %[none] shift t to put it in third dimestion
    t1 =permute(t,[3 2 1]);
    erf=(1i./pi).*trapz(t1,exp(-t.^2)./(z-t),3);
end