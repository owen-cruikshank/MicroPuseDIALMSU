% Calculating pertabative absorption of O2
% author: Owen Cruikshank
% date: 11/30/2020

function [alpha_final,alpha_1_raw,alpha_2_raw,Spectrum] = pertAbsorption(alpha, T_etalon, Model, Range, Time, Spectrum, BSR, ind_r_lo,ind_r_hi, Options, UseRBspectrum,usePCAabsorption)
    % --- Spectral distribution using the initial temperature profile guess ---

     cB = 1.2;%Brullouion correction to doppler gaussian half width
%     cB = -0.01*(Range.rm+altitude) + 1.2;%Brullouin correction for 1.2 at 0km and 1.1 at 10km
     cB = -0.01*((Range.rm+1500)/1000) + 1.2;%Brullouin correction for 1.2 at 0km and 1.1 at 10km
% 
%     c_doppler_O2 = m_air*c^2./(8*(nu_online*100).^2*kb);                   %[m^2 K] Doppler coefficient
%     doppler_O2_un_ret = ((c_doppler_O2./T/pi).^0.5).*exp(-c_doppler_O2.*(nu_online*100-nu_scan_3D_short*100).^2./T./cB.^2); %[m] Doppler broadended lineshape         
% 
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
    dopplerWidth = dopplerWidth.*29.9792458;


    norm_O2_ret = trapz(doppler_O2_un_ret,3).*Spectrum.nuBin*100;                   %[none] Lineshape integral
    doppler_O2_ret = doppler_O2_un_ret./norm_O2_ret;                       %[m] Normalized doppler lineshape


        c_doppler_O2 = m_air*c^2./(8*(Spectrum.nu_offline(1)*100).^2*kb);                   %[m^2 K] Doppler coefficient
    doppler_O2_un_ret = ((c_doppler_O2./Model.T/pi).^0.5).*exp(-c_doppler_O2.*(Spectrum.nu_offline(1)*100-Spectrum.nu_scan_3D_short_off*100).^2./Model.T./cB.^2); %[m] Doppler broadended lineshape         

    norm_O2_ret = trapz(doppler_O2_un_ret,3).*Spectrum.nuBin*100;                   %[none] Lineshape integral
    doppler_O2_ret_off = doppler_O2_un_ret./norm_O2_ret;                       %[m] Normalized doppler lineshape

    clear norm_O2_ret c_doppler_O2 doppler_O2_un_ret
    % Check if doppler_o2_ret is normalized to 1 when integrated across frequency
    %doppler_o2_ret_check = trapz(doppler_O2_ret,3).*nuBin*100;              %[none]


% % % 
% % % T_err = 3;
% % % 
% % % T_err = 3;

% cl = sqrt((m_air*c^2)/(8*(Spectrum.nu_online(1)*100).^2*kb)/pi);
% c2 = -(m_air*c^2)/(8*(Spectrum.nu_online(1)*100).^2*kb) .*(Spectrum.nu_online(1)*100-Spectrum.nu_scan_3D_short*100).^2 ./cB.^2;
% dldt= cl.*(-0.5).*Model.T.^(-3/2).*exp(c2./Model.T)+cl./Model.T.*c2./Model.T.*(-Model.T.^2).*exp(c2./Model.T);
% l_err = sqrt((T_err.*dldt).^2);%error setimate for l


% % % l_err = T_err .* abs(sqrt(c_doppler_O2./pi).*exp(-c_doppler_O2.*(Spectrum.nu_online(1)*100-Spectrum.nu_scan_3D_short*100).^2./Model.T./cB.^2) .*(-0.5.*Model.T.^(-3/2)+Model.T.^(-0.5).*(    c_doppler_O2.*(Spectrum.nu_online(1)*100-Spectrum.nu_scan_3D_short*100).^2./Model.T.^2./cB.^2)));
% % % l_err = l_err./norm_O2_ret;
% % % clear c1 c2 dldt
% % % 
% % % deltafunc = zeros(size(Spectrum.nu_scan_3D_short));
% % % deltafunc(:,:,Spectrum.online_index(1)) = 1; %create delta function for aerosol scattering
% % % 


    %%
    %Calculate RB spectrum by PCA
else
    [doppler_O2_ret] = RB_O2_770_PCA(Model.T,Model.P,Spectrum.nu_scan_3D_short,Spectrum);  
end
    %%
  
% % %     BSR_err = 0.10 * abs(HSRL.BSR);
% % %    % BSR_err = .03 * abs(HSRL.BSR);
% % % %%g1_err = sqrt((BSR_err.*(-HSRL.BSR.^(-2).*doppler_O2_ret + HSRL.BSR.^(-2).*deltafunc)).^2 + (l_err./HSRL.BSR).^2); %error estimate for g
% % % 
% % % g1_err = sqrt((BSR_err.*HSRL.BSR.^(-2).*(deltafunc-doppler_O2_ret)).^2 + (l_err.*HSRL.BSR.^(-1)).^2);
% % % 
% % % clear l_err deltafunc

    % --- Backscatter Lineshape g ---
    g1_m = 1./BSR .* doppler_O2_ret ;%.*nuBin*100;                         %[m] Molecular backscatter lineshape
    %%clear doppler_O2_ret
    g1_a = zeros(Range.i_range,Time.i_time,Spectrum.i_scan_3D_short);                       % Initalize aerosol lineshape
    for i = 1:Time.i_time
        g1_a(:,i,Spectrum.online_index(1)) = (1 - 1./BSR(:,i))/ Spectrum.nuBin / 100 ; %[m] aerosol backscatter lineshape
    end
    g1 = g1_a + g1_m;                                                   %[m] Combined backscatter lineshape
    %g1_check = trapz(g1,3).*nuBin*100;                                %[none] Check if integral of g1 is normalized to 1  

%     g1_m = 1./BSR .* doppler_O2_ret_off ;%.*nuBin*100; 
%     g2 = g1_a+g1_m;
    
    clear g1_m g1_a
    %derivative of lineshape dg/dr
    ind_r_lo = 1:(length(Range.rm)-1);
    ind_r_hi = 2:length(Range.rm);
    dg1_dr = (g1(ind_r_hi,:,:) - g1(ind_r_lo,:,:)) ./(Range.rangeBin*Options.oversample); %[none] Derivative over oversamped range
    %dg1_dr = interp1(Range.rm(ind_r_lo),dg1_dr,Range.rm,'nearest',nan);         %[none] Make dg/dr the same size as g
   dg1_dr(ind_r_hi(end),:,:) = dg1_dr(ind_r_hi(end-1),:,:);

%     for iii = 1:length(Spectrum.lambda_scan_3D_short)
%         dg1_dr1(:,:,iii) = interp2(Time.ts,Range.rm(ind_r_lo)+Range.rangeBin./2,dg1_dr(:,:,iii),Time.ts,Range.rm);
%         dg1_dr1(:,:,iii) = fillmissing(dg1_dr1(:,:,iii),'nearest');
%     end
%     dg1_dr = dg1_dr1;
% 
%     %first order central derivative
    dg1 = zeros(size(g1));
    dg1(1,:,:) = (g1(2,:,:)-g1(1,:,:))./Range.rangeBin;
    for iii = 2:length(Range.rm)-1
        dg1(iii,:,:) = (g1(iii+1,:,:)-g1(iii-1,:,:))/2/Range.rangeBin;
    end
    dg1(end,:,:) = (g1(end,:,:)-g1(end-1,:,:))./Range.rangeBin;

    dg1_dr = dg1;

%     dg2 = zeros(size(g2));
%     dg2(1,:,:) = (g2(2,:,:)-g2(1,:,:))./Range.rangeBin;
%     for iii = 2:length(Range.rm)-1
%         dg2(iii,:,:) = (g2(iii+1,:,:)-g2(iii-1,:,:))/2/Range.rangeBin;
%     end
%     dg2(end,:,:) = (g2(end,:,:)-g2(end-1,:,:))./Range.rangeBin;
% 
%     dg2_dr = dg2;

    clear dg1


% % %     for iii = 1:size(dg1_dr,1)-1
% % %         dg1_dr_err(iii,:,:) = sqrt(g1_err(iii,:,:).^2+g1_err(iii+1,:,:).^2)./(Range.rangeBin*Options.oversample);
% % %     end
% % %     dg1_dr_err(iii+1,:,:)=dg1_dr_err(iii,:,:); %derivative error
    %%
%       ==old absorption function
%     %[absorption_f] = absorption_O2_770_model_wavenumber(T,P,nu_scan_3D_short,WV); %[m] lineshape function 
%     for i = 1:i_range
%         parfor j=1:i_time
%     absorption_f(i,j,:) = absorption_O2_770_model_wavenumber(T(i,j),P(i,j),nu_scan_3D_short(1,1,:),WV(i,j)); %[m-1] Funcrtion to calculate theoretical absorption
%         end
%     end

    %%
    disp('PCA absorption')
  %  if UseRBspectrum
   % [absorption_f,~] = absorption_O2_770_PCA(Model.T,Model.P,Spectrum.nu_scan_3D_short(1,1,:),Model.WV);
   % 
   %  else
   %  absorption_f = absorption_O2_770_model_wavenumber(Model.T,Model.P,Spectrum.nu_scan_3D_short ,Model.WV);
   %  end

   absorption_f = absorption_O2_770_PCA(Model.T,Model.P,Spectrum.nu_scan_3D_short ,Model.WV);



    %o2absorption_off = absorption_O2_770_model_wavenumber(Model.T,Model.P,Spectrum.nu_scan_3D_short_off ,Model.WV);


    o2absorption_off = zeros(size(absorption_f));
    o2absorption_off = absorption_O2_770_model(Model.T,Model.P,Spectrum.nu_offline,Model.WV);
    %%
    %Create lineshape function
    for i = 1:Time.i_time
        absorption_f(:,i,:) = absorption_f(:,i,:) ./ absorption_f(:,i,Spectrum.online_index(1));  %[none] Normalize lineshape function
       % absorption_f(:,i,:) = absorption_f(:,i,:) ./ max(absorption_f(:,i,:),3);  %[none] Normalize lineshape function
    end

% % %     absorption_f_err = abs(0.02*absorption_f); %error for absorption model
   %%% alpha_err = abs(0.10*alpha);

     %%    
    % --- Zeroth Order Transmission ---
    Tm0 = exp(-cumtrapz(Range.rm,alpha.*absorption_f,1));      %[none] Zeroth order transmission  
    %Tm0 = exp(-cumtrapz(Range.rm,absorption_f,1));      %[none] Zeroth order transmission  


    TmOff = exp(-cumtrapz(Range.rm,o2absorption_off,1)); 
% % %     atimesAF_err = abs(alpha.*absorption_f).*sqrt((alpha_err./alpha).^2+(absorption_f_err./absorption_f).^2);


%     cumtrap_err = atimesAF_err(1,:,:).^2;
%     for iii = 2:size(atimesAF_err,1)
%         cumtrap_err = atimesAF_err(iii,:,:).^2 + cumtrap_err;
%     end

% % % for iii = 1:size(atimesAF_err,2)
% % %     for jjj = 1:size(atimesAF_err,3)
% % %         for kkk = 1:size(atimesAF_err,1)
% % %              cumtrap_err(kkk,iii,jjj)=sumsqr(atimesAF_err(1:kkk,iii,jjj));
% % %         end
% % %     end
% % % end

    
% % %     cumtrap_err = sqrt(cumtrap_err).*Range.rangeBin;
% % %     clear atimesAF_err
% % % 
% % %     Tm0_err = sqrt((cumtrap_err.*(-1).*Tm0).^2);
% % %     clear cumtrap_err

    % Integrand terms
    % Online
    zeta = g1.*T_etalon;                        %[m]
    eta = dg1_dr.*T_etalon;                     %[none]
% % % 
% % %     zeta_err = g1_err.*T_etalon;
% % %     eta_err = dg1_dr_err.*T_etalon;

    clear g1 dg1_dr
    % Integrated terms
    % Online
    zeta_int = trapz(zeta.*Tm0,3);              %[none]
    %eta_int = trapz(eta.*Tm0,3)*Spectrum.nuBin*100;                %[1/m]
    %zeta_ls_int = trapz(zeta.*Tm0.*(1-absorption_f),3)*Spectrum.nuBin*100;    %[none]
    % Offline
    %zeta2_int = trapz(zeta,3)*Spectrum.nuBin*100;                  %[none]
    %eta2_int = trapz(eta,3)*Spectrum.nuBin*100;                    %[1/m]

    % === First Order ===
   % W1 = zeta_ls_int./zeta_int;                 %[none]
    %W1 = trapz(zeta.*Tm0.*(1-absorption_f),3)./trapz(zeta.*Tm0,3);                 %[none]
    W1 = trapz(zeta.*Tm0.*(1-absorption_f),3)./zeta_int;    
    
% % %     W1insidenum_err = abs(zeta.*Tm0.*(1-absorption_f)) .* sqrt((zeta_err./zeta).^2+(Tm0_err./Tm0).^2+(absorption_f_err./(1-absorption_f)).^2);
% % %     W1insideden_err = abs(zeta.*Tm0) .* sqrt((zeta_err./zeta).^2+(Tm0_err./Tm0).^2);
% % % 
% % %     G1insidenum_err = abs(eta.*Tm0) .* sqrt((eta_err./eta).^2+(Tm0_err./Tm0).^2);

%     W1_err_num = abs(zeta.*Tm0.*(1-absorption_f)) .* sqrt((g1_err.*T_etalon./zeta).^2 + (Tm0_err./Tm0).^2 + (absorption_f_err./(1-absorption_f)).^2);
% 
%     W1_err_denom = abs(zeta.*Tm0) .* sqrt((g1_err.*T_etalon./zeta).^2 + (Tm0_err./Tm0).^2);
% 
%     G1_err_num = abs(eta.*Tm0) .* sqrt((dg1_dr_err.*T_etalon./eta).^2 + (Tm0_err./Tm0).^2);
% 
%     G1_err_num1 = G1_err_num(:,:,1).^2;
%     W1_err_num1 = W1_err_num(:,:,1).^2;
%     W1_err_denom1 = W1_err_denom(:,:,1).^2;
% % % W1insidenum_err1  = W1insidenum_err(:,:,1).^2;
% % % W1insideden_err1 = W1insideden_err(:,:,1).^2;
% % % G1insidenum_err1 = G1insidenum_err(:,:,1).^2;
% % %     for iii = 2:size(Spectrum.nu_scan_3D_short)
% % % %         W1_err_num1 = W1_err_num(:,:,iii).^2+W1_err_num1;
% % % %         W1_err_denom1 = W1_err_denom(:,:,iii).^2+W1_err_denom1;
% % % %         G1_err_num1 = G1_err_num(:,:,iii).^2+G1_err_num1;
% % % W1insidenum_err1 = W1insidenum_err(:,:,iii).^2 +W1insidenum_err1;
% % % W1insideden_err1 = W1insideden_err(:,:,iii).^2 +W1insideden_err1;
% % % G1insidenum_err1 =G1insidenum_err(:,:,iii).^2 +G1insidenum_err1;
% % %     end
% % % 
% % % W1insidenum_err1 = sqrt(W1insidenum_err1);
% % % W1insideden_err1 =sqrt(W1insideden_err1);
% % % G1insidenum_err1 =sqrt(G1insidenum_err1);
% % % 
% % % W1_err = abs(W1).*sqrt((W1insidenum_err1./(trapz(zeta.*Tm0.*(1-absorption_f),3))).^2+(W1insideden_err1./(trapz(zeta.*Tm0,3))).^2);
% % % 
% % % G1_err = abs(trapz(eta.*Tm0,3)./trapz(zeta.*Tm0,3)).*sqrt((G1insidenum_err1./(trapz(eta.*Tm0,3))).^2+(W1insideden_err1./(trapz(zeta.*Tm0,3))).^2);
    
%     clear W1_err_denom W1_err_num G1_err_num
%     W1_err_num1 = sqrt(W1_err_num1);
%     W1_err_denom1 = sqrt(W1_err_denom1);
%     G1_err_num1 = sqrt(G1_err_num1);
% 
% 
%     W1_err = abs(W1) .* sqrt((W1_err_num1./trapz(zeta.*Tm0.*(1-absorption_f),3)).^2 + (W1_err_denom1./trapz(zeta.*Tm0,3)).^2);
%     clear W1_err_num1 
% 
%     G1_err = abs(trapz(eta.*Tm0,3)./trapz(zeta.*Tm0,3)) .* sqrt((G1_err_num1./trapz(eta.*Tm0,3)).^2 + (W1_err_denom1./trapz(zeta.*Tm0,3)).^2);
%     clear G1_err_num1 W1_err_denom1
     

    %G1 = eta_int./zeta_int - eta2_int./zeta2_int;%[1/m]
    %G1 = eta_int./zeta_int - trapz(zeta,3)./trapz(eta,3);%[1/m]
    %G1 = trapz(eta.*Tm0,3)./trapz(zeta.*Tm0,3) - trapz(eta,3)./trapz(zeta,3);%[1/m]

    %alpha_1_raw = 0.5.*(alpha.*W1 + G1);      %[1/m]

   % alpha_1_raw = 0.5.*(alpha.*W1 +  trapz(eta.*Tm0,3)./trapz(zeta.*Tm0,3) - trapz(eta,3)./trapz(zeta,3));      %[1/m]
    alpha_1_raw = 0.5.*(alpha.*W1 +  trapz(eta.*Tm0,3)./zeta_int - trapz(eta.*TmOff,3)./trapz(zeta.*TmOff,3));      %[1/m]

%     alpha_1_raw2 = 0.5.*(alpha.*W1 +  trapz(eta.*Tm0,3)./zeta_int - trapz(dg2_dr.*T_etalon.*TmOff,3)./trapz(g2.*T_etalon.*TmOff,3));      %[1/m]

% % %     alpha_1_raw_err = sqrt((W1_err.*0.5.*alpha).^2+(alpha_err.*0.5.*W1).^2 + (G1_err.*0.5).^2);
    % --- First Order Transmission Tm1 ---

% nanAlpha = isnan(alpha_1_raw);
% alpha_1_raw(nanAlpha) = Model.absorption;

alpha_1_raw=fillmissing(alpha_1_raw,'nearest',1);

    Tm1 = exp(-cumtrapz(Range.rm,Options.oversample.*alpha_1_raw.*absorption_f,1));      %[none] First order transmission

    % === Second Order ===
    % Integrated terms
    zeta_Tm1_int = trapz(zeta.*Tm0.*(1-Tm1),3)*Spectrum.nuBin*100;             %[none]
    eta_Tm1_int = trapz(eta.*Tm0.*(1-Tm1),3)*Spectrum.nuBin*100;               %[1/m]
    zeta_ls_Tm1_int = trapz(zeta.*Tm0.*(1-Tm1).*(1-absorption_f),3)*Spectrum.nuBin*100;   %[none]
    
    clear Tm1

%     W2 = (zeta_ls_int.*zeta_ls_Tm1_int./(zeta_int.^2)) - (zeta_ls_Tm1_int./zeta_int);   %[none]
%     G2 = (eta_int.*zeta_Tm1_int./(zeta_int.^2)) - (eta_Tm1_int./zeta_int);              %[1/m]

   % W2 = ((trapz(zeta.*Tm0.*(1-absorption_f),3)*Spectrum.nuBin*100).*zeta_ls_Tm1_int./((trapz(zeta.*Tm0,3)*Spectrum.nuBin*100).^2)) - (zeta_ls_Tm1_int./(trapz(zeta.*Tm0,3)*Spectrum.nuBin*100));   %[none]
    W2 = ((trapz(zeta.*Tm0.*(1-absorption_f),3)*Spectrum.nuBin*100).*zeta_ls_Tm1_int./((zeta_int*Spectrum.nuBin*100).^2)) - (zeta_ls_Tm1_int./(zeta_int*Spectrum.nuBin*100));   %[none]
 
    %G2 = ((trapz(eta.*Tm0,3)*Spectrum.nuBin*100).*zeta_Tm1_int./((trapz(zeta.*Tm0,3)*Spectrum.nuBin*100).^2)) - (eta_Tm1_int./(trapz(zeta.*Tm0,3)*Spectrum.nuBin*100));              %[1/m]
    G2 = ((trapz(eta.*Tm0,3)*Spectrum.nuBin*100).*zeta_Tm1_int./((zeta_int*Spectrum.nuBin*100).^2)) - (eta_Tm1_int./(zeta_int*Spectrum.nuBin*100));              %[1/m]

    alpha_2_raw = 0.5.*(alpha_1_raw.*W1 + alpha.*W2 + G2);    %[1/m]
    
    alpha_final = alpha_2_raw+alpha_1_raw+alpha;
    

%     Spectrum.g = doppler_O2_ret;
%     Spectrum.g1 = g1;
%     Spectrum.l = absorption_f;

end