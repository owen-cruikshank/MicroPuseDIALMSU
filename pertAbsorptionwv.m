% Calculating pertabative absorption of O2
% author: Owen Cruikshank
% date: 2/14/2022

% function [alpha_1, alpha_2,Spectrum] = pertAbsorptionwv(alpha_0, T_etalon, T, P, rm,ts,rkm,m_air, nu_online, nu_scanwv_3D_short,nuBin,BSR,ind_r_lo,ind_r_hi,WV,online_index,i_range,i_time,i_scan_3D_short,rangeBin,oversample,t_avg,c,kb,altitude,Spectrum)
%     function [alpha_1, alpha_2,Spectrum] = pertAbsorptionwv(alpha_0, T_etalon,Model, Range, Time,Constant,Spectrum, nu_online, nu_scanwv_3D_short,nuBin,BSR,ind_r_lo,ind_r_hi,WV,online_index,i_range,i_time,i_scan_3D_short,rangeBin,oversample,t_avg,c,kb,altitude,Spectrum)
 function [alpha_1, alpha_2,Spectrum] = pertAbsorptionwv(alpha_0, T_etalon, Model, Range, Time, Spectrum, BSR828, ind_r_lo,ind_r_hi, Options, Constant, altitude)
    % --- Spectral distribution using the initial temperature profile guess ---

    cB = -0.01*(Range.rkm+altitude) + 1.2;%Brullouin correction for 1.2 at 0km and 1.1 at 10km

    c_doppler_O2 = Constant.m_air*Constant.c^2./(8*(Spectrum.nu_wvon*100).^2*Constant.kb);                   %[m^2 K] Doppler coefficient
    doppler_O2_ret = ((c_doppler_O2./Model.T/pi).^0.5).*exp(-c_doppler_O2.*(Spectrum.nu_wvon*100-Spectrum.nu_scanwv_3D_short*100).^2./Model.T./cB.^2); %[m] Doppler broadended lineshape         

    clear c_doppler_O2

    norm_O2_ret = trapz(doppler_O2_ret,3).*Spectrum.nuBin*100;                   %[none] Lineshape integral
    doppler_O2_ret = doppler_O2_ret./norm_O2_ret;                       %[m] Normalized doppler lineshape

    clear norm_O2_ret
    % Check if doppler_o2_ret is normalized to 1 when integrated across frequency
    %doppler_o2_ret_check = trapz(doppler_O2_ret,3).*nuBin*100;              %[none]


    %%
    %Calculate RB spectrum by PCA
    %[doppler_O2_ret] = RB_O2_770_PCA(T,P,nu_scan_3D_short);  
% % %     [doppler_O2_ret] = RB_828_PCA(Model.T,Model.P,Spectrum.nu_scanwv_3D_short);
% % % 
% % %     %shift to correct online wavelength
% % %     RBshift = 151-Spectrum.online_indexwv;
% % %     RBshift = 145-Spectrum.online_indexwv;
% % %     RBshift = zeros(size(Spectrum.online_indexwv));
% % %     RBshift = 152-Spectrum.online_indexwv;
% % %     for iii=1:length(Spectrum.online_indexwv)
% % %         if RBshift(iii)<=0
% % %             doppler_O2_ret(:,iii,:) = cat(3,zeros(Range.i_range,1,-RBshift(iii)),doppler_O2_ret(:,iii,1:end+RBshift(iii)));
% % %         elseif RBshift(iii)>0
% % %             doppler_O2_ret(:,iii,:) = cat(3,doppler_O2_ret(:,iii,RBshift(iii)+1:end),zeros(Range.i_range,1,RBshift(iii)));
% % %         end
% % %     end

     norm_O2_ret = trapz(doppler_O2_ret,3).*Spectrum.nuBin*100;                   %[none] Lineshape integral
    doppler_O2_ret = doppler_O2_ret./norm_O2_ret;                       %[m] Normalized doppler lineshape

    %%
  
    % --- Backscatter Lineshape g ---
    g1_m = 1./BSR828 .* doppler_O2_ret ;%.*nuBin*100;                         %[m] Molecular backscatter lineshape
    g1_a = zeros(Range.i_range,Time.i_time,size(Spectrum.nu_scanwv_3D_short,3));                       % Initalize aerosol lineshape
    for i = 1:Time.i_time
        g1_a(:,i,Spectrum.online_indexwv(i)) = (1 - 1./BSR828(:,i))/ Spectrum.nuBin / 100 ; %[m] aerosol backscatter lineshape
    end
    g1 = g1_a + g1_m;                                                   %[m] Combined backscatter lineshape
    %g1_check = trapz(g1,3).*nuBin*100;                                %[none] Check if integral of g1 is normalized to 1

    clear g1_a g1_m doppler_O2_ret
    %derivative of lineshape dg/dr
    dg1_dr = (g1(ind_r_hi,:,:) - g1(ind_r_lo,:,:)) ./(Range.rangeBin*Options.oversample); %[none] Derivative over oversamped range
    dg1_dr = interp1(Range.rm(ind_r_lo),dg1_dr,Range.rm,'nearest',nan);         %[none] Make dg/dr the same size as g

    %%
%       ==old absorption function

% disp('Calculation absorption linewidth')
%     absorption_f = zeros(length(T(:,1)),length(T(1,:)),length(nu_scanwv_3D_short(1,1,:)));
%     %for i = 1:i_range
%         for j=1:i_time
%             %absorption_f(:,j,:) = cross_section_wv_828_model_wavenumber(T(:,j),P(:,j),nu_scanwv_3D_short(1,1,:));
%         end
%     %end

    %%
     disp('PCA absorption')
    absorption_f = cross_section_wv_828_PCA(Model.T,Model.P,Spectrum.nu_scanwv_3D_short(1,1,:));

    %%
    %Create lineshape function
    f = nan(size(absorption_f));
    for i = 1:Time.i_time
        f(:,i,:) = absorption_f(:,i,:) ./ absorption_f(:,i,Spectrum.online_indexwv(i));  %[none] Normalize lineshape function
    end

    clear absorption_f
     %%    
    % --- Zeroth Order Transmission ---
    Tm0 = exp(-cumtrapz(Range.rm,alpha_0.*f,1));      %[none] Zeroth order transmission  

    % Integrand terms
    % Online
    zeta = g1.*T_etalon;                        %[m]
    eta = dg1_dr.*T_etalon;                     %[none]


    clear dg1_dr g1
    % Integrated terms
    % Online
    zeta_int = trapz(zeta.*Tm0,3)*Spectrum.nuBin*100;              %[none]
    eta_int = trapz(eta.*Tm0,3)*Spectrum.nuBin*100;                %[1/m]
    zeta_ls_int = trapz(zeta.*Tm0.*(1-f),3)*Spectrum.nuBin*100;    %[none]
    % Offline
    zeta2_int = trapz(zeta,3)*Spectrum.nuBin*100;                  %[none]
    eta2_int = trapz(eta,3)*Spectrum.nuBin*100;                    %[1/m]

    % === First Order ===
    W1 = zeta_ls_int./zeta_int;                 %[none]
    G1 = eta_int./zeta_int - eta2_int./zeta2_int;%[1/m]

    alpha_1 = 0.5.*(alpha_0.*W1 + G1);      %[1/m]


    % --- First Order Transmission Tm1 ---
    Tm1 = exp(-cumtrapz(Range.rm,Options.oversample.*alpha_1.*f,1));      %[none] First order transmission

    % === Second Order ===
    % Integrated terms
    zeta_Tm1_int = trapz(zeta.*Tm0.*(1-Tm1),3)*Spectrum.nuBin*100;             %[none]
    eta_Tm1_int = trapz(eta.*Tm0.*(1-Tm1),3)*Spectrum.nuBin*100;               %[1/m]
    zeta_ls_Tm1_int = trapz(zeta.*Tm0.*(1-Tm1).*(1-f),3)*Spectrum.nuBin*100;   %[none]

    W2 = (zeta_ls_int.*zeta_ls_Tm1_int./(zeta_int.^2)) - (zeta_ls_Tm1_int./zeta_int);   %[none]
    G2 = (eta_int.*zeta_Tm1_int./(zeta_int.^2)) - (eta_Tm1_int./zeta_int);              %[1/m]

    alpha_2 = 0.5.*(alpha_1.*W1 + alpha_0.*W2 + G2);    %[1/m]


%     Spectrum.gwv = doppler_O2_ret;
%     Spectrum.g1wv = g1;
%     Spectrum.lwv = absorption_f;

end