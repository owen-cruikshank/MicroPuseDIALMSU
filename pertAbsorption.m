% Calculating pertabative absorption of O2
% author: Owen Cruikshank
% date: 11/30/2020

function [alpha_final,alpha_1_raw,alpha_2_raw,Spectrum] = pertAbsorption(alpha, T_etalon, T_etalon_off, Model, Range, Time, Spectrum, BSR, Options, UseRBspectrum,usePCAabsorption,Constants)
    % --- Spectral distribution using the initial temperature profile guess ---

     cB = 1.2;%Brullouion correction to doppler gaussian half width
%     cB = -0.01*(Range.rm+altitude) + 1.2;%Brullouin correction for 1.2 at 0km and 1.1 at 10km
     cB = -0.01*((Range.rm+1500)/1000) + 1.2;%Brullouin correction for 1.2 at 0km and 1.1 at 10km
% 
%     c_doppler_O2 = m_air*c^2./(8*(nu_online*100).^2*kb);                   %[m^2 K] Doppler coefficient
%     doppler_O2_un_ret = ((c_doppler_O2./T/pi).^0.5).*exp(-c_doppler_O2.*(nu_online*100-nu_scan_3D_short*100).^2./T./cB.^2); %[m] Doppler broadended lineshape         
%     norm_O2_ret = trapz(doppler_O2_un_ret,3).*nuBin*100;                   %[none] Lineshape integral
%     doppler_O2_ret = doppler_O2_un_ret./norm_O2_ret;                       %[m] Normalized doppler lineshape
% 
%     % Check if doppler_o2_ret is normalized to 1 when integrated across frequency
%     %doppler_o2_ret_check = trapz(doppler_O2_ret,3).*nuBin*100;              %[none]
% 
if ~UseRBspectrum
    m_air = 28.97/1000./6.02214e23;
    kb = 1.3806e-23;
    c=3e8;

    c_doppler_O2 = m_air*c^2./(8*(Spectrum.nu_online(1)*100).^2*kb);                   %[m^2 K] Doppler coefficient
    doppler_O2_un_ret = ((c_doppler_O2./Model.T/pi).^0.5).*exp(-c_doppler_O2.*(Spectrum.nu_online(1)*100-Spectrum.nu_scan_3D_short*100).^2./Model.T./cB.^2); %[m] Doppler broadended lineshape         
    dopplerWidth = sqrt(Model.T./c_doppler_O2./log(2)).*cB;%FWHM cm^-1]

    norm_O2_ret = trapz(doppler_O2_un_ret,3).*Spectrum.nuBin*100;                   %[none] Lineshape integral
    doppler_O2_ret = doppler_O2_un_ret./norm_O2_ret;                       %[m] Normalized doppler lineshape

    c_doppler_O2 = m_air*c^2./(8*(Spectrum.nu_offline(1)*100).^2*kb);                   %[m^2 K] Doppler coefficient
    doppler_O2_un_ret = ((c_doppler_O2./Model.T/pi).^0.5).*exp(-c_doppler_O2.*(Spectrum.nu_offline(1)*100-Spectrum.nu_scan_3D_short_off*100).^2./Model.T./cB.^2); %[m] Doppler broadended lineshape         

    norm_O2_ret = trapz(doppler_O2_un_ret,3).*Spectrum.nuBin*100;                   %[none] Lineshape integral
    doppler_O2_ret_off = doppler_O2_un_ret./norm_O2_ret;                       %[m] Normalized doppler lineshape

    clear norm_O2_ret c_doppler_O2 doppler_O2_un_ret
    % Check if doppler_o2_ret is normalized to 1 when integrated across frequency
    %doppler_o2_ret_check = trapz(doppler_O2_ret,3).*nuBin*100;              %[none]
    
else
    %Calculate RB spectrum by PCA
    [doppler_O2_ret_offline] = RB_O2_770_PCA(Model.T,Model.P,Spectrum.nu_scan_3D_short_off,Spectrum);  
    [doppler_O2_ret_online] = RB_O2_770_PCA_online(Model.T,Model.P,Spectrum.nu_scan_3D_short,Spectrum); 
end
    %%

    % --- Backscatter Lineshape g ---
    g1on_m = 1./BSR .* doppler_O2_ret_online ;%.*nuBin*100;                         %[m] Molecular backscatter lineshape
    g1on_a = zeros(Range.i_range,Time.i_time,Spectrum.i_scan_3D_short); % Initalize aerosol lineshape
    g1off_m = 1./BSR .* doppler_O2_ret_offline ;%.*nuBin*100;                         %[m] Molecular backscatter lineshape
    g1off_a = zeros(Range.i_range,Time.i_time,Spectrum.i_scan_3D_short);  
    
    for i = 1:Time.i_time
        g1on_a(:,i,Spectrum.online_index(1)) = (1 - 1./BSR(:,i))/ Spectrum.nuBin / 100 ; %[m] aerosol backscatter lineshape
        g1off_a(:,i,Spectrum.offline_index(1)) = (1 - 1./BSR(:,i))/ Spectrum.nuBin / 100 ; %[m] aerosol backscatter lineshape
    end
    g1on = g1on_a + g1on_m;                                                   %[m] Combined backscatter lineshape
    g1off = g1off_a + g1off_m; 
    %g1_check = trapz(g1,3).*nuBin*100;                                %[none] Check if integral of g1 is normalized to 1  
    
    clear g1on_m g1on_a g1off_m g1off_a
    %derivative of lineshape dg/dr
   %  ind_r_lo = 1:(length(Range.rm)-1);
   %  ind_r_hi = 2:length(Range.rm);
   %  dg1_dr = (g1(ind_r_hi,:,:) - g1(ind_r_lo,:,:)) ./(Range.rangeBin*Options.oversample); %[none] Derivative over oversamped range
   %  %dg1_dr = interp1(Range.rm(ind_r_lo),dg1_dr,Range.rm,'nearest',nan);         %[none] Make dg/dr the same size as g
   % dg1_dr(ind_r_hi(end),:,:) = dg1_dr(ind_r_hi(end-1),:,:);

%     for iii = 1:length(Spectrum.lambda_scan_3D_short)
%         dg1_dr1(:,:,iii) = interp2(Time.ts,Range.rm(ind_r_lo)+Range.rangeBin./2,dg1_dr(:,:,iii),Time.ts,Range.rm);
%         dg1_dr1(:,:,iii) = fillmissing(dg1_dr1(:,:,iii),'nearest');
%     end
%     dg1_dr = dg1_dr1;
% 
    %first order central derivative
    dg1on = zeros(size(g1on));
    dg1off = zeros(size(g1off));
    dg1on(1,:,:) = (g1on(2,:,:)-g1on(1,:,:))./Range.rangeBin;
    dg1off(1,:,:) = (g1off(2,:,:)-g1off(1,:,:))./Range.rangeBin;
    for iii = 2:length(Range.rm)-1
        dg1on(iii,:,:) = (g1on(iii+1,:,:)-g1on(iii-1,:,:))/2/Range.rangeBin;
        dg1off(iii,:,:) = (g1off(iii+1,:,:)-g1off(iii-1,:,:))/2/Range.rangeBin;
    end
    dg1on(end,:,:) = (g1on(end,:,:)-g1on(end-1,:,:))./Range.rangeBin;
    dg1off(end,:,:) = (g1off(end,:,:)-g1off(end-1,:,:))./Range.rangeBin;


    %%
    disp('PCA absorption')

    %absorption_f = absorption_O2_770_PCA(Model.T,Model.P,Spectrum.nu_scan_3D_short ,Model.WV);
    absorption_f = absorption_O2_770_PCA2(Model.T,Model.P,Spectrum,Model.WV,Constants);
    absorption = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_online,Model.WV,Constants);

    o2absorption_off = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_offline,Model.WV,Constants);
    %%
    %Create lineshape function
    for i = 1:Time.i_time
        %absorption_f(:,i,:) = absorption_f(:,i,:) ./ absorption_f(:,i,Spectrum.online_index(1));  %[none] Normalize lineshape function
        absorption_f(:,i,:) = absorption_f(:,i,:) ./ absorption(:,i);  %[none] Normalize lineshape function
    end

     %%    
    % --- Zeroth Order Transmission ---
    Tm0 = exp(-cumtrapz(Range.rm,alpha.*absorption_f,1));      %[none] Zeroth order transmission  
    TmOff = exp(-cumtrapz(Range.rm,o2absorption_off,1)); 

    % Integrand terms
    % Online
    zeta = g1on.*T_etalon;                        %[m]
    eta = dg1on.*T_etalon;                     %[none]

    zeta_off = g1off.*T_etalon_off;
    eta_off = dg1off.*T_etalon_off;

    % zeta_off = g1off.*T_etalon;
    % eta_off = dg1off.*T_etalon;

    clear g1on g1off dg1on dg1off
    % Integrated terms

    % === First Order ===  
    W1 = trapz(zeta.*Tm0.*(1-absorption_f),3)./trapz(zeta.*Tm0,3); 
    

    alpha_1_raw = 0.5.*(alpha.*W1 +  trapz(eta.*Tm0,3)./trapz(zeta.*Tm0,3) - trapz(eta_off.*TmOff,3)./trapz(zeta_off.*TmOff,3));      %[1/m]

    % --- First Order Transmission Tm1 ---
    alpha_1_raw=fillmissing(alpha_1_raw,'nearest',1);

    Tm1 = exp(-cumtrapz(Range.rm,Options.oversample.*alpha_1_raw.*absorption_f,1));      %[none] First order transmission

    % === Second Order ===
    % Integrated terms
    zeta_Tm1_int = trapz(zeta.*Tm0.*(1-Tm1),3)*Spectrum.nuBin*100;             %[none]
    eta_Tm1_int = trapz(eta.*Tm0.*(1-Tm1),3)*Spectrum.nuBin*100;               %[1/m]
    zeta_ls_Tm1_int = trapz(zeta.*Tm0.*(1-Tm1).*(1-absorption_f),3)*Spectrum.nuBin*100;   %[none]
    
    clear Tm1

    W2 = ((trapz(zeta.*Tm0.*(1-absorption_f),3)*Spectrum.nuBin*100).*zeta_ls_Tm1_int./((trapz(zeta.*Tm0,3)*Spectrum.nuBin*100).^2)) - (zeta_ls_Tm1_int./(trapz(zeta.*Tm0,3)*Spectrum.nuBin*100));   %[none]
 
    G2 = ((trapz(eta.*Tm0,3)*Spectrum.nuBin*100).*zeta_Tm1_int./((trapz(zeta.*Tm0,3)*Spectrum.nuBin*100).^2)) - (eta_Tm1_int./(trapz(zeta.*Tm0,3)*Spectrum.nuBin*100));              %[1/m]

    alpha_2_raw = 0.5.*(alpha_1_raw.*W1 + alpha.*W2 + G2);    %[1/m]
    %== total alpha
    alpha_final = alpha_2_raw+alpha_1_raw+alpha;

%     Spectrum.g = doppler_O2_ret;
%     Spectrum.g1 = g1;
%     Spectrum.l = absorption_f;
end