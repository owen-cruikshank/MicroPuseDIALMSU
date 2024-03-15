function plot_O2(p_point,sonde_index,~,Sonde,Model,Counts,Range,Time,~,Temperature,Format,Alpha,~,HSRL,Data,~,cloud_SDm_above,~,~,N_wvm,N_wv0m,~,~)

%function plot_O2(p_point,sonde_index,span_days,Sonde,Model,Counts,Range,Time,~,Temperature,Format,Alpha,cloud_SDm,HSRL,Data,SNRm,cloud_SDm_above,N_wv,N_wv0,N_wvm,N_wv0m,AbsHumm,AbsHum0m)

%==Fitted lapse rate and surface temperature
figure(728590)
subplot(2,1,1)
plot(Time.date_ts,permute(Temperature.L_fit_sm_test(1,:,end)*1000,[3 2 1]))
%hold on
%plot(Time.ts([p_point(1) p_point(end)])/60/60,[-10 -4],'k','linewidth',2)
%plot(Time.date_ts,movmean(permute(Temperature.L_fit_sm_test(1,:,end)*1000,[3 2 1]),40),'linewidth',2)
%hold off
yline(-6.5);
grid on
title('Fitted lapse rate')
xlabel('Time (UTC hr)')
ylabel('Lapse rate (K/km)')
%legend('Raw','Smoothed')
subplot(2,1,2)
%plot(Time.date_ts,-Model.Ts(1,:) + permute(Temperature.Ts_fit(1,:,:),[3 2 1]))
%hold on
plot(Time.date_ts,-Model.Ts(1,:) + permute(Temperature.Ts_fit(1,:,end),[3 2 1]))
%%plot(Time.date_ts,movmean(-Model.Ts(1,:) + permute(Temperature.Ts_fit(1,:,end),[3 2 1]),40))
%plot(Time.ts([p_point(1) p_point(end)])/60/60,[-20 20])
%hold off
grid on
title('Surface temp - fitted surface temp')
xlabel('Time (UTC hr)')
ylabel('\Delta(T_{fit}-T_{surf})')


figure(8983745)
if ~isempty(Sonde.sonde_ind)
subplot(2,1,1)
plot(Time.date_ts(1,Sonde.sonde_ind(1,:)),Temperature.L_fit_sm_test(1,Sonde.sonde_ind(1,:),end)*1000,'-*')
grid on
ylabel('Lapse (K/km)')
title('Fitted Lapse rate at radiosonde launches')
subplot(2,1,2)
plot(Time.date_ts(1,Sonde.sonde_ind(1,:)),-Model.Ts(1,Sonde.sonde_ind(1,:)) + Temperature.Ts_fit(1,Sonde.sonde_ind(1,:),end),'-*')
grid on
ylabel('\Delta(T_{fit}-T_{surf})')
end




%==Background Counts
figure(9985)
%semilogy(datenum(Time.date_ts),Counts.bg_o2off,'-',datenum(Time.date_ts),Counts.bg_o2on,datenum(Time.date_ts),Counts.bg_o2off_mol,datenum(Time.date_ts),Counts.bg_o2on_mol,datenum(Time.date_ts),Counts.bg_wvon,datenum(Time.date_ts),Counts.bg_wvoff,datenum(Time.date_ts([p_point(1) p_point(end)])),[1 100])
semilogy(Time.date_ts,Counts.bg_o2off,'--')
    hold on
   semilogy(Time.date_ts,Counts.bg_o2on,'--')
    semilogy(Time.date_ts,Counts.bg_o2off_mol)
    semilogy(Time.date_ts,Counts.bg_o2on_mol)
    semilogy(Time.date_ts,Counts.bg_wvon)
    semilogy(Time.date_ts,Counts.bg_wvoff)
    
hold on
for i = 1:size(Sonde.sonde_ind,2)
    semilogy([Time.date_ts(:,Sonde.sonde_ind(1,i)) Time.date_ts(:,Sonde.sonde_ind(end,i))],[1 10^3],'k')
end
hold off
%timeAxis(Format,Time)
title('Background')
xlabel('Time UTC hr')
ylabel('photons')
legend('Offline','Online','Offline molecular','Online molecular','Online WV','Offline WV')
grid on

%=Diff Counts and sonde difference
% figure(99496)
% cla reset
% line(diff(Counts.o2off(:,p_point(1)))/Range.rangeBin,Range.rkm(1:end-1),'Color','r')
% line(diff(Counts.o2off(:,p_point(1)))./Range.rangeBin./Counts.o2off(1:end-1,p_point(1))*100,Range.rkm(1:end-1),'Color','b')
% %line(diag(dBSRSG(:,p_point)./BSR(:,p_point)),Range.rkm,'Color','b')
% %line(diag(dBSRlo(:,p_point)),Range.rkm(1:end-1))
% legend('dOFF/dr')
% xlim([-2 2])
% xlabel('dOFF/dr')
% ylabel('Range (km)')
% ax1 = gca; % current axes
% ax1.XColor = 'r';
% ax1.YColor = 'r';
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% %line(diag(alpha_0(:,p_point)-Sonde.absorption_sonde{sonde_index}),Range.rkm,'Parent',ax2,'Color','k')
% line(diag(Alpha.alpha_0_filt(:,p_point)-Sonde.absorption_sonde{sonde_index}),Range.rkm,'Parent',ax2,'Color','g')
% %line((diag(alpha_0_filt(1:end-1,p_point)))./(diff(Counts.o2off(:,p_point(1)))./Range.rangeBin./Counts.o2off(1:end-1,p_point(1))*100)-Sonde.absorption_sonde{sonde_index}(1:end-1),Range.rkm(1:end-1),'Parent',ax2,'Color','k')
% %line((diag(alpha_0_filt(1:end-1,p_point)))./(diff(Counts.o2off(:,p_point(1)))./Range.rangeBin)-Sonde.absorption_sonde{sonde_index}(1:end-1),Range.rkm(1:end-1),'Parent',ax2,'Color','k')
% grid on
% xlim([-.5 .5]*10^-3)
% xlabel('\alpha_0 - \alpha_{sonde}')
% ylabel('Range (km)')
% legend('\alpha_0 - \alpha_{sonde}')


%=Water Vapor profile
figure(99410)
plot(diag(N_wv0m(:,p_point)),Range.rkm,'Linewidth',2)
hold on
plot(diag(N_wvm(:,p_point)),Range.rkm,'--','Linewidth',2)
%plot(diag(N_wv0(:,p_point)),Range.rkm)
plot(Sonde.WV_sonde(:,sonde_index),Range.rkm,'.-')
plot(diag(Model.WV(:,p_point)),Range.rkm)
hold off
xlim([1e22 3e23])
ylim([0 4])
grid on
legend('VW0','VW puturbative','Radiosonde')
xlabel('WV (molec/m^3)')
ylabel('Range (km)')

%=Temperature profile
figure(884)
%plot(Ts_fit(1,p_point,end),0,'+')
%temp = diag(Temperature.T_final_test(:,p_point))-273.13;
% plot(temp,Range.rkm)
%plot(diag(Temperature.T_final_tests0(:,p_point))-273.13,Range.rkm,'linewidth',2)
%hold on
plot(diag(Temperature.T_finalm(:,p_point))-273.13,Range.rkm,'linewidth',2)
hold on
%plot(temp,Range.rkm)
%Temperature.exclusion(Temperature.exclusion==0)=NaN;

plot(Sonde.T_sonde(:,sonde_index)-273.13,Range.rm/1000,'.-')

% plot(diag(Temperature.T_finalm(:,p_point))-273.13+diag(Temperature.TempStd(:,p_point)),Range.rkm,'--')
% plot(diag(Temperature.T_finalm(:,p_point))-273.13-diag(Temperature.TempStd(:,p_point)),Range.rkm,'--')

%%plot(diag(Temperature.T_finalm(:,p_point))-273.13+diag(Temperature.TempStds(:,p_point)),Range.rkm,'*-')
%%plot(diag(Temperature.T_finalm(:,p_point))-273.13-diag(Temperature.TempStds(:,p_point)),Range.rkm,'*-')
 plot(diag(Model.T(:,p_point))-273.13,Range.rkm,'-','LineWidth',2)
%%plot(Sonde.T_sonde(:,sonde_index)+2,Range.rm/1000,'--k')
%%plot(Sonde.T_sonde(:,sonde_index)-2,Range.rm/1000,'--k')
plot(diag(Temperature.L_fit_sm_test(:,p_point,end) .* Range.rm + Temperature.Ts_fit(1,p_point,end))-273.13,Range.rkm,'--','linewidth',2)

plot(Sonde.T_sonde(:,sonde_index)-273.13+2,Range.rm/1000,'.-')
plot(Sonde.T_sonde(:,sonde_index)-273.13-2,Range.rm/1000,'.-')
plot(diag(Temperature.T_final_test(:,p_point))-273.13,Range.rkm,'--')
hold off
grid on
xlabel('Temperature (^oC)')
ylabel('Range (km)')
ylim([0 4])
xlim([-15 20])
title(sprintf(['Temperature\n' datestr(Time.date_ts(p_point(1)))]))
legend('Retrieved temperature','Radiosonde','model','fit')

%=Tsonde-T
figure(888)
plot(Sonde.T_sonde(:,sonde_index) - diag(Temperature.T_finalm(:,p_point)),Range.rkm)
%plot(diag(Model.T(:,p_point)) - diag(T_finalm(:,p_point)),Range.rkm)
hold on
plot(Sonde.T_sonde(:,sonde_index) - diag(Model.T(:,p_point)),Range.rkm)

plot(Sonde.T_sonde(:,sonde_index) - diag(Temperature.T_finalm(:,p_point)-Temperature.TempStds(:,p_point)),Range.rkm,'--')
plot(Sonde.T_sonde(:,sonde_index) - diag(Temperature.T_finalm(:,p_point)+Temperature.TempStds(:,p_point)),Range.rkm,'--')
plot(diag(Temperature.TempStds(:,p_point)),Range.rkm,'*-')
%plot(diag(Temperature.L_fit_sm_test(:,p_point,end) .* Range.rm + Temperature.Ts_fit(1,p_point,end))-diag(Temperature.T_finalm(:,p_point)),Range.rkm)
hold off
line([0 0],[0 6],'Color','k','LineStyle','--')
line([-1 -1],[0 6],'Color','k','LineStyle','-.')
line([1 1],[0 6],'Color','k','LineStyle','-.','HandleVisibility','off')
line([-2 -2],[0 6],'Color','k')
line([2 2],[0 6],'Color','k')
xlabel('\Delta T (T_{sonde} - T_{retrieved}) (K)')
ylabel('Range (km)')
ylim([0 4])
xlim([-8 8])
title(sprintf(['Temperature difference\n' datestr(Time.date_ts(p_point(1)))]))
%legend('MPD','lapse','fit')

%===WV percent difference
% % figure(888999)
% % plot((Sonde.WV_sonde(:,sonde_index) - diag(N_wv0(:,p_point)))./Sonde.WV_sonde(:,sonde_index)*100,Range.rkm)
% % %plot(diag(Model.T(:,p_point)) - diag(T_finalm(:,p_point)),Range.rkm)
% % hold on
% % plot((Sonde.WV_sonde(:,sonde_index) - diag(N_wvm(:,p_point)))./Sonde.WV_sonde(:,sonde_index)*100,Range.rkm)
% % plot((Sonde.WV_sonde(:,sonde_index) - diag(N_wv0m(:,p_point)))./Sonde.WV_sonde(:,sonde_index)*100,Range.rkm)
% % hold off
% % line([0 0],[0 6],'Color','k','LineStyle','--')
% % line([-10 -10],[0 6],'Color','k','LineStyle','-.')
% % line([10 10],[0 6],'Color','k','LineStyle','-.','HandleVisibility','off')
% % line([-20 -20],[0 6],'Color','k')
% % line([20 20],[0 6],'Color','k')
% % xlabel('(WV_{sonde} - WV_{retrieved})/WV_{sonde}*100 (%)')
% % ylabel('Range (km)')
% % ylim([0 5])
% % xlim([-100 100])
% % title(sprintf(['WVdifference\n' datestr(Time.date_ts(p_point(1)))]))
% % legend('MPD 0','MPD total','MPD0')

%=2D temperature
% % figure(886)
% % r_max_plot = 4; % Max range to plot [km]
% % [~,ind_km_max] = min(abs(Range.rkm-r_max_plot)); % Max range Index
% % imAlpha=ones(size(Temperature.T_finalm(1:ind_km_max,:))); % initalize transpernecy matrix
% % imAlpha(isnan(Temperature.T_finalm(1:ind_km_max,:)))=0; % Set AlphaData to not plot NaNs
% % imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1; %Set transperency so that cloud mask can be seen
% % colorLimits=[-25 20];
% % imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),Temperature.T_finalm(1:ind_km_max,:)-273.13,'AlphaData',imAlpha,colorLimits) %Plot
% % hold on
% % for ii = 1:size(Sonde.sonde_ind,2)
% %      plot(datenum(Time.date_ts([Sonde.sonde_ind(1,ii) Sonde.sonde_ind(end,ii)])),[Range.rkm(1) Range.rkm(end)],'--k')
% % end
% % plot(datenum([Time.date_ts(p_point(1)) Time.date_ts(p_point(1))]),[Range.rkm(1) Range.rkm(end)],'--k')
% % colors = colormap; % Set colors to colormap
% % colors(1,:) = [0 0 0]; % Set lowest color black to make clouds appear black
% % colormap(colors); % Make colormap new colors
% % set(gca, 'YDir','normal') % Set y axis increasing
% % set(gca,'color',[1 1 1]);% Color background white
% % colorbar % Add colorbar
% % timeAxis(Format,Time)
% % title(sprintf('Temperature retrieval')) %Add title
% % ylabel('Range (km)')
% % title(colorbar,'Temperature (^oC)') % Add title to colorbar
% % hold off

% % % 
% % % figure(8866666)
% % % r_max_plot = 5; % Max range to plot [km]
% % % [~,ind_km_max] = min(abs(Range.rkm-r_max_plot)); % Max range Index
% % % surf(Time.date_ts,Range.rkm(1:ind_km_max),Temperature.T_finalm(1:ind_km_max,:)-273,'linestyle','none')
% % % colorbar % Add colorbar
% % % view(2)
% % % title(sprintf('Temperature retrieval')) %Add title
% % % ylabel('Range (km)')
% % % title(colorbar,'Temperature (^oC)') % Add title to colorbar
% % % hold off


%===2D temp with error mask
figure(888666)
h = pcolor(Time.date_ts,Range.rkm,Temperature.T_finalm-273);
set(h, 'EdgeColor', 'none');
hold on
for ii = 1:size(Sonde.sonde_ind,2)
    plot(Time.date_ts(Sonde.sonde_ind(:,ii)),Range.rkm,'--k')
end
colorbar % Add colorbar
clim([-5 35])
%colormap("hot")
colormap("turbo")
title(sprintf('Temperature retrieval')) %Add title
ylabel('Range (km)')
title(colorbar,'Temperature (^oC)') % Add title to colorbar
ylim([ 0 4])
hold off

figure(887)
%h = pcolor(Time.date_ts,Range.rkm,Model.WV);
h = pcolor(Time.date_ts,Range.rkm,Alpha.AbsHum0Rawm);
set(h, 'EdgeColor', 'none');
hold on
for ii = 1:size(Sonde.sonde_ind,2)
    plot(Time.date_ts(Sonde.sonde_ind(:,ii)),Range.rkm,'--k')
end
colorbar % Add colorbar
clim([0 3e23])
clim([0 12])
%colormap("hot")
title(sprintf('WV')) %Add title
ylabel('Range (km)')
title(colorbar,'WV (g/m^3)') % Add title to colorbar
ylim([ 0 4])
hold off

%===2D water vapor
% % figure(887)
% % r_max_plot = 4; % Max range to plot [km]
% % [~,ind_km_max] = min(abs(Range.rkm-r_max_plot)); % Max range Index
% % % imAlpha=ones(size(Temperature.T_finalm(1:ind_km_max,:))); % initalize transpernecy matrix
% % % imAlpha(isnan(Temperature.T_finalm(1:ind_km_max,:)))=0; % Set AlphaData to not plot NaNs
% % % imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=0; %Set transperency so that cloud mask can be seen
% % %imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1; %Set transperency so that cloud mask can be seen
% % colorLimits =[0 3e23];
% % imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),N_wvm(1:ind_km_max,:)-273.13,colorLimits) %Plot
% % colormap('jet')
% % colors = colormap; % Set colors to colormap
% % colors(1,:) = [0 0 0]; % Set lowest color black to make clouds appear black
% % colormap(colors); % Make colormap new colors
% % set(gca, 'YDir','normal') % Set y axis increasing
% % set(gca,'color',[1 1 1]);% Color background white
% % colorbar % Add colorbar
% % timeAxis(Format,Time)
% % title(sprintf('WV number density')) %Add title
% % xlabel('Time')
% % ylabel('Range (km)')
% % title(colorbar,'WV (1/m^3)') % Add title to colorbar

% figure(886645)
% r_max_plot = 4; % Max range to plot [km]
% [~,ind_km_max] = min(abs(Range.rkm-r_max_plot)); % Max range Index
% colorLimits =[-1 1]*10^-3;
% imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),movstd(Alpha.alpha_total_rawf(1:ind_km_max,:)-Alpha.alpha_total_rawg(1:ind_km_max,:),5,0,1),colorLimits) %Plot
% hold on
% %sonde lines
% % %  plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b','LineWidth',2)
% % % %yline(.6)
% % % for ii = 1:size(Sonde.sonde_ind,2)
% % %     plot(datenum(Time.date_ts([Sonde.sonde_ind(1,ii) Sonde.sonde_ind(end,ii)])),[Range.rkm(1) Range.rkm(end)],'--k')
% % % end
% colormap('jet')
% colors = colormap; % Set colors to colormap
% colors(1,:) = [0 0 0]; % Set lowest color black to make clouds appear black
% colormap(colors); % Make colormap new colors
% colormap(redblue)
% set(gca, 'YDir','normal') % Set y axis increasing
% %set(gca,'color',[1 1 1]);% Color background white
% colorbar % Add colorbar
% timeAxis(Format,Time)
% title(sprintf('Temperature retrieval')) %Add title
% ylabel('Range (km)')
% title(colorbar,'Temperature (^oC)') % Add title to colorbar
% hold off

%=Absorption profile
figure(4)
%alpha = Alpha.alpha_0;
alpha = Alpha.alpha_0_full;
%alpha(cloud_SDm_above)=nan;
%alpha(mask)=nan;

alpha_tot = Alpha.alpha_totalm;
alpha_tot(cloud_SDm_above)=nan;
%alpha_tot(mask)=nan;
%plot(diag(alpha(:,p_point)),Range.rkm)
%plot(diag(Alpha.alpha_0s(:,p_point)),Range.rkm,'LineWidth',2)
plot(diag(alpha(:,p_point)),Range.rkm,'LineWidth',2)
hold on
plot(diag(alpha_tot(:,p_point)),Range.rkm,'LineWidth',2)
plot(Sonde.absorption_sonde{sonde_index},Range.rkm,'.-')
% plot(diag(Alpha.alpha_totalm(:,p_point)+sqrt(Alpha.alpha_total_err(:,p_point).^2)),Range.rkm,'.-','LineWidth',1)
% plot(diag(Alpha.alpha_totalm(:,p_point)-sqrt(Alpha.alpha_total_err(:,p_point).^2)),Range.rkm,'.-','LineWidth',1)
% plot(diag(Alpha.alpha_totalm(:,p_point)+sqrt(Alpha.alpha_total_errs(:,p_point).^2)),Range.rkm,'*-','LineWidth',1)
% plot(diag(Alpha.alpha_totalm(:,p_point)-sqrt(Alpha.alpha_total_errs(:,p_point).^2)),Range.rkm,'*-','LineWidth',1)
% plot(diag(Alpha.alpha_totalm(:,p_point)+(Alpha.alpha_0_err(:,p_point))),Range.rkm,'--','LineWidth',1)
% plot(diag(Alpha.alpha_totalm(:,p_point)-(Alpha.alpha_0_err(:,p_point))),Range.rkm,'--','LineWidth',1)
%plot(diag(Alpha.alpha_total_raw(:,p_point)),Range.rkm,'--','LineWidth',2)
%plot(diag(Alpha.alpha_0m(:,p_point)),Range.rkm,'LineWidth',2)
%plot(diag(Model.absorption(:,p_point)),Range.rkm)
%plot(diag(alpha_d(:,p_point)),Range.rkm,'--','LineWidth',2)
%plot(diag(alpha_0_filt(:,p_point)),Range.rkm,'--','LineWidth',2)
%plot(diag(alpha_0_filt_mol(:,p_point)),Range.rkm,'--','LineWidth',2)
%plot(alpha_0_filt_decon(:,p_point(1)),Range.rkm,'.-')
%plot(alpha_d(:,p_point(1)),Range.rkm(1:end-1),'.-')
%plot(alpha_dModel(:,p_point(1)),Range.rkm(1:end-1),'.-')
hold off
grid on
title(sprintf(['O_2 absorption coefficient\n' datestr(Time.date_ts(p_point(1))) 'UTC']))
ylabel('Range (km)')
xlabel('Absorption (m^{-1})')
%xlim([-2e-4 6e-4])
xlim([0 2e-4])
ylim([0 4])
legend('Zero order','Smoothed corrected','Radiosonde model')

%=WV Absorption profile
% % figure(5)
% % plot(diag(Alpha.alpha_0wv(:,p_point)),Range.rkm);
% % hold on
% % plot(diag(Alpha.alpha_total_rawwv(:,p_point)),Range.rkm);
% % hold on
% % hold off
% % grid on
% % title(sprintf(['Calculated WV absorption coefficient\n' datestr(Time.date_ts(p_point(1)))]))
% % ylabel('Range (km)')
% % xlabel('Absorption (m^{-1})')
% % xlim([-2e-4 6e-4])
% % ylim([0 5])
% % legend('Measured zeroth order absorption','Purtabative absorption')


%=BSR
figure(444)
bsr = HSRL.BSR;
bsr(cloud_SDm_above)=nan;
plot(diag(bsr(:,p_point)),Range.rkm,'linewidth',2)
hold on
%plot(diag(HSRL.BSR(:,p_point)),Range.rkm)
%%%plot(diag(HSRL.BSR828(:,p_point)),Range.rkm,'linewidth',2)
hold off
grid on
title(sprintf(['Backscatter ratio\n' datestr(Time.date_ts(p_point(1)))]))
ylabel('Range (km)')
xlim([0 3])
ylim([0 4])

%figure(445)
%  dg1 = zeros(size(HSRL.BSR));
%     dg1(1,:,:) = (HSRL.BSR(2,:,:)-HSRL.BSR(1,:,:))./Range.rangeBin;
%     for iii = 2:length(Range.rm)-1
%         dg1(iii,:,:) = (HSRL.BSR(iii+1,:,:)-HSRL.BSR(iii-1,:,:))/2/Range.rangeBin;
%     end
%     dg1(end,:,:) = (HSRL.BSR(end,:,:)-HSRL.BSR(end-1,:,:))./Range.rangeBin;
%     
%     plot(diag(-dg1(:,p_point)).*1000,Range.rkm)
%     hold off


%=Counts Profile
figure(7)
subplot(2,2,1)
plot(diag(Counts.o2on(:,p_point)),Range.rm)
hold on
plot(diag(Counts.o2off(:,p_point)),Range.rm)
plot(diag(Counts.o2on_mol(:,p_point)),Range.rm)
plot(diag(Counts.o2off_mol(:,p_point)),Range.rm)
plot(Counts.o2on_bgsub(:,p_point(1)),Range.rm_raw_o2,'--')
plot(Counts.o2off_bgsub(:,p_point(1)),Range.rm_raw_o2,'--')
plot(Counts.o2on_bgsub_mol(:,p_point(1)),Range.rm_raw_o2,'--')
plot(Counts.o2off_bgsub_mol(:,p_point(1)),Range.rm_raw_o2,'--')
xlim([0 500])
title('Averaged photon returns')
legend('Online','Offline','Online molecular','Offline molecular','raw online','raw offline','raw online molecular','raw offline molecular')
hold off
grid on
ylabel('Range [m]')
xlabel('Photons')
subplot(2,2,2)
plot(log(Counts.o2on_bgsub(:,p_point(1)).*Range.rm_raw_o2.^2),Range.rm_raw_o2/1000);
hold on
plot(log(Counts.o2off_bgsub(:,p_point(1)).*Range.rm_raw_o2.^2),Range.rm_raw_o2/1000);
plot(log(Counts.o2on_bgsub_mol(:,p_point(1)).*Range.rm_raw_o2.^2),Range.rm_raw_o2/1000);
plot(log(Counts.o2off_bgsub_mol(:,p_point(1)).*Range.rm_raw_o2.^2),Range.rm_raw_o2/1000);
legend('on','off','on mol','off mol')
xlabel('log(Photons*range^2)')
ylabel('Range [m]')
grid on
hold off
subplot(2,2,3)
semilogx((Counts.o2on_bgsub(:,p_point(1))),Range.rm_raw_o2/1000);
hold on
semilogx((Counts.o2off_bgsub(:,p_point(1))),Range.rm_raw_o2/1000);
semilogx((Counts.o2on_bgsub_mol(:,p_point(1))),Range.rm_raw_o2/1000);
semilogx(Counts.o2off_bgsub_mol(:,p_point(1)),Range.rm_raw_o2/1000);
legend('on','off','on mol','off mol')
xlabel('Photons')
ylabel('Range [m]')
grid on
hold off
subplot(2,2,4)
semilogx(Counts.o2on_bgsub(:,p_point(1))+Counts.bg_o2on(:,p_point(1)),Range.rm_raw_o2/1000);
hold on
semilogx(Counts.o2off_bgsub(:,p_point(1))+Counts.bg_o2off(:,p_point(1)),Range.rm_raw_o2/1000);
semilogx(Counts.o2on_bgsub_mol(:,p_point(1))+Counts.bg_o2on_mol(:,p_point(1)),Range.rm_raw_o2/1000);
semilogx(Counts.o2off_bgsub_mol(:,p_point(1))+Counts.bg_o2off_mol(:,p_point(1)),Range.rm_raw_o2/1000);
legend('on','off','on mol','off mol')
xlabel('Photons')
ylabel('Range [m]')
grid on
hold off

%=Housekeeping temperature
% % figure(738732)
% % subplot(4,1,1)
% % plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.InsideCell.Temperature)
% % legend('Inside Cell')
% % title('Thermocouples')
% % ylabel('^o C')
% % subplot(4,1,2)
% % plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.OutsideCell.Temperature  )
% % ylabel('^o C')
% % legend('Outside Cell')
% % subplot(4,1,3)
% % plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.TSOA.Temperature  )
% % ylabel('^o C')
% % legend('TSOA')
% % subplot(4,1,4)
% % plot(Data.Thermocouple.InsideCell.TimeStamp  ,Data.Thermocouple.RoomTemp.Temperature  )
% % ylabel('^o C')
% % legend('Room Temp')
% % xlabel('Time UTC hr')

%=Laser locking wavelenths
figure(738734)
subplot(3,1,1)
plot(Data.Laser.O2Online.TimeStamp,Data.Laser.O2Online.WaveDiff)
hold on
plot(Data.Laser.O2Online.TimeStamp,Data.Laser.O2Offline.WaveDiff)
hold off
legend('Online','Offline')
title('wavelength Difference')
ylim([-5*10^-4 5*10^-4])
ylabel('(nm)')
xlabel('Time hr')
subplot(3,1,2)
plot(Data.Laser.O2Online.TimeStamp  ,Data.Laser.O2Online.WaveDiff-Data.Laser.O2Offline.WaveDiff  )
title('Online wave diff-offline wave diff')
ylim([-5*10^-4 5*10^-4])
grid on
ylabel('(nm)')
xlabel('Time hr')
subplot(3,1,3)
%=Laser seed power
%figure(778298)
plot(Data.Laser.O2Online.TimeStamp,10.^(Data.Laser.O2Online.SeedPower/10))
hold on
plot(Data.Laser.O2Offline.TimeStamp,10.^(Data.Laser.O2Offline.SeedPower/10))
hold off
title('Laser Seed power at wavemeter')
legend('Online','Offline')

%=Laser temperatures
% % figure(738735)
% % subplot(4,1,1)
% % plot(Data.Etalon.O2Etalon.TimeStamp,Data.Laser.O2Online.TemperatureActual)
% % hold on
% % plot(Data.Etalon.O2Etalon.TimeStamp,Data.Laser.O2Online.TemperatureDesired)
% % legend('Actual','Desired')
% % hold off
% % grid on
% % title('Online Laser Temperature')
% % ylabel('(C)')
% % subplot(4,1,2)
% % plot(Data.Etalon.O2Etalon.TimeStamp,Data.Laser.O2Offline.TemperatureActual)
% % hold on
% % plot(Data.Etalon.O2Etalon.TimeStamp,Data.Laser.O2Offline.TemperatureDesired)
% % legend('Actual','Desired')
% % hold off
% % grid on
% % title('Offline Laser Temperature')
% % ylabel('(C)')
% % subplot(4,1,3)
% % % plot(Data.Laser.TWSOA.TimeStamp,Data.Laser.TWSOA.TemperatureActual)
% % % hold on
% % % grid on
% % % plot(Data.Laser.TWSOA.TimeStamp,Data.Laser.TWSOA.TemperatureDesired)
% % % legend('Actual','Desired')
% % % hold off
% % % title('TWSOA Laser Temperature')
% % % ylabel('(C)')
% % subplot(4,1,4)
% % plot(Data.Etalon.O2Etalon.TimeStamp,Data.Etalon.O2Etalon.TemperatureActual)
% % title('Etalon Temperature')
% % grid on
% % ylabel('(C)')
% % xlabel('Time hr')



%=Weather Station
figure(83)
subplot(2,1,1)
plot(Time.date_ts,Model.Ts-273.13)
hold on
plot(Time.date_ts,Temperature.Ts_fit(1,:,end)-273.13)
for i=1:size(Sonde.sonde_ind,2)
plot(Time.date_ts(Sonde.sonde_ind(1,i)),Sonde.Tsurf(1,i)-273.13,'*')
end
hold off
legend('Weather station','Fit')
%timeAxis(Format,Time)
ylabel('^oC')
%xlim(Time.date_ts(1)) datenum(Time.date_ts(end))])
subplot(2,1,2)
plot(Time.date_ts,Model.Ps)
%timeAxis(Format,Time)
title('Weather Station Suface Pressure')
xlabel('Time')
ylabel('Pressure ATM')
%xlim([datenum(Time.date_ts(1)) datenum(Time.date_ts(end))])

% % %=UPS
% % if ishandle(648)
% %     figure(648)
% %     clf(648)
% % end
% % figure(648)
% % line(Data.UPS.all.TimeStamp,Data.UPS.all.BatteryTimeLeft,'Color','r')
% % line(Data.UPS.all.TimeStamp,Data.UPS.all.HoursOnBattery,'Color','g')
% % ax1=gca;
% % ax1.XColor = 'r';
% % ax1.YColor = 'r';
% % xlabel('time (hrs)')
% % ylabel('Battery remaining (hrs)')
% % legend('Battery remaining','Hours on battery','Location','Southwest')
% % ax1_pos = ax1.Position; % position of first axes
% % ax2 = axes('Position',ax1_pos,...
% %     'XAxisLocation','top',...
% %     'YAxisLocation','right',...
% %     'Color','none');
% % line(Data.UPS.all.TimeStamp,Data.UPS.all.UPSTemperature,'Parent',ax2,'Color','k')
% % ylabel('Battery temperature (C)')

%=2d delta
% figure(4489)
% imagesc(Time.thr,Range.rm/1000,Counts.delta)
% hold on
% xline(Time.thr(p_point(1)),'K');
% hold off
% set(gca, 'YDir','normal')
% caxis([-1000 1000]/3)
% colorbar
% colormap(redblue(64))
% title('Counts.delta')
% ylabel('Range (km)')
% xlabel('Time (UTC hrs)')

%=2D counts
% % figure(589)
% % subplot(2,2,1)
% % imagesc(Time.thr,Range.rm/1000,Counts.o2on_noise.*Range.rm.^2)
% % hold on
% % xline(Time.thr(p_point(1)),'r');
% % hold off
% % set(gca, 'YDir','normal')
% % caxis([0 1e9])
% % colorbar
% % colormap('hot')
% % title('o2 on unaveraged * r^{ 2}')
% % ylabel('Range (km)')
% % xlabel('Time (UTC hrs)')
% % subplot(2,2,2)
% % imagesc(Time.thr,Range.rm/1000,Counts.o2on_noise)
% % hold on
% % xline(Time.thr(p_point(1)),'r');
% % hold off
% % set(gca, 'YDir','normal')
% % caxis([0 500])
% % colorbar
% % colormap('hot')
% % title('o2 on unaveraged ')
% % ylabel('Range (km)')
% % xlabel('Time (UTC hrs)')
% % subplot(2,2,3)
% % imagesc(Time.thr,Range.rm,Counts.o2off_noise.*Range.rm.^2)
% % hold on
% % xline(Time.thr(p_point(1)),'r');
% % hold off
% % set(gca, 'YDir','normal')
% % caxis([0 1e9])
% % colorbar
% % title('o2off unaveraged  * r^{ 2}')
% % subplot(2,2,4)
% % imagesc(Time.thr,Range.rm,Counts.o2off_noise)
% % hold on
% % xline(Time.thr(p_point(1)),'r');
% % hold off
% % set(gca, 'YDir','normal')
% % caxis([0 500])
% % colorbar
% % title('o2off unaveraged')


%=2D molecular counts
figure(590)
subplot(2,2,1)
imagesc(Time.thr,Range.rm/1000,Counts.o2on_noise_mol.*Range.rm.^2)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
clim([0 2.5e11])
colorbar
colormap('hot')
title('o2 on mol unaveraged * r^{ 2}')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,2)
imagesc(Time.thr,Range.rm/1000,Counts.o2on_noise_mol)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
clim([0 2e5])
colorbar
colormap('hot')
title('o2 on mol unaveraged')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,3)
imagesc(Time.thr,Range.rm,Counts.o2off_noise_mol.*Range.rm.^2)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
clim([0 2.5e11])
colorbar
colormap('hot')
title('o2 off mol unaveraged * r^{ 2}')
subplot(2,2,4)
imagesc(Time.thr,Range.rm,Counts.o2off_noise_mol)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
clim([0 0.5e5])
colorbar
colormap('hot')
title('o2 off mol unaveraged')

figure(591)
subplot(2,2,1)
imagesc(Time.thr,Range.rm/1000,Counts.o2on_noise.*Range.rm.^2)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
clim([0 2.5e11])
colorbar
colormap('hot')
title('o2 on unaveraged * r^{ 2}')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,2)
imagesc(Time.thr,Range.rm/1000,Counts.o2on_noise)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
clim([0 2e5])
colorbar
colormap('hot')
title('o2 on unaveraged')
ylabel('Range (km)')
xlabel('Time (UTC hrs)')
subplot(2,2,3)
imagesc(Time.thr,Range.rm,Counts.o2off_noise.*Range.rm.^2)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
clim([0 2.5e11])
colorbar
colormap('hot')
title('o2 off  unaveraged * r^{ 2}')
subplot(2,2,4)
imagesc(Time.thr,Range.rm,Counts.o2off_noise)
hold on
xline(Time.thr(p_point(1)),'r');
hold off
set(gca, 'YDir','normal')
clim([0 0.5e5])
colorbar
colormap('hot')
title('o2 off unaveraged')

%===unaveraged counts
% figure(590)
% subplot(2,1,1)
% imagesc(Time.thr,Range.rm/1000,Counts.o2on_noise_mol)
% hold on
% xline(Time.thr(p_point(1)),'r');
% hold off
% set(gca, 'YDir','normal')
% set(gca,'ColorScale','log')
% caxis([0 1e3])
% colorbar
% colormap('hot')
% title('o2 on mol unaveraged')
% ylabel('Range (km)')
% xlabel('Time (UTC hrs)')
% subplot(2,1,2)
% imagesc(Time.thr,Range.rm,Counts.o2off_noise_mol)
% hold on
% xline(Time.thr(p_point(1)),'r');
% hold off
% set(gca, 'YDir','normal')
% set(gca,'ColorScale','log')
% caxis([0 100])
% colorbar
% colormap('hot')
% title('o2 off mol unaveraged')

%===WV counts
% figure(591)
% subplot(2,1,1)
% imagesc(Time.thr,Range.rm/1000,Counts.wvon_noise)
% hold on
% xline(Time.thr(p_point(1)),'r');
% hold off
% set(gca, 'YDir','normal')
% set(gca,'ColorScale','log')
% caxis([0 1e3])
% colorbar
% colormap('hot')
% title('wv on unaveraged')
% ylabel('Range (km)')
% xlabel('Time (UTC hrs)')
% subplot(2,1,2)
% imagesc(Time.thr,Range.rm,Counts.wvoff_noise)
% hold on
% xline(Time.thr(p_point(1)),'r');
% hold off
% set(gca, 'YDir','normal')
% set(gca,'ColorScale','log')
% caxis([0 1000])
% colorbar
% colormap('hot')
% title('wv off unaveraged')

%=2D absorption
figure(8886)
subplot(4,1,1)
r_max_plot = 6; % Max range to plot [km]
[~,ind_km_max] = min(abs(Range.rkm-r_max_plot));
colorLimits = [-1e-4, 3e-4];
h=pcolor(Time.date_ts,Range.rkm(1:ind_km_max),Alpha.alpha_0(1:ind_km_max,:));
set(h, 'EdgeColor', 'none');
hold on
clim(colorLimits)
colorbar
title(sprintf('O2 alpha 0 %s to %s (UTC)',Time.date_ts(1),Time.date_ts(end)))
title('Zero order DIAL O_2 absorption')
ylabel('Range (km)')
title(colorbar,'m^{-1}')
set(gca, 'XTickLabel', [])
hold off
subplot(4,1,2)
colorLimits = [-5e-5, 5e-5];
h=pcolor(Time.date_ts,Range.rkm(1:ind_km_max),Alpha.alpha_1(1:ind_km_max,:));
set(h, 'EdgeColor', 'none');
hold on
clim(colorLimits)
colorbar
title(sprintf('alpha 1'))
ylabel('Range (km)')
title(colorbar,'m^{-1}')
set(gca, 'XTickLabel', [])
hold off
subplot(4,1,3)
colorLimits = [-3e-6, 3e-6];
h=pcolor(Time.date_ts,Range.rkm(1:ind_km_max),Alpha.alpha_2(1:ind_km_max,:));
set(h, 'EdgeColor', 'none');
hold on
clim(colorLimits)
colorbar
title(sprintf('alpha 2'))
ylabel('Range (km)')
title(colorbar,'m^{-1}')
set(gca, 'XTickLabel', [])
hold off
subplot(4,1,4)
colorLimits = [-1e-4, 3e-4];
h=pcolor(Time.date_ts,Range.rkm(1:ind_km_max),Alpha.alpha_total_raw(1:ind_km_max,:));
set(h, 'EdgeColor', 'none');
hold on
clim(colorLimits)
colors = colormap;
colors(1,:) = [0 0 0];%set lowest color black
colormap(colors);
colorbar
title('alpha total')
ylabel('Range (km)')
title(colorbar,'m^{-1}')
xlabel('Time UTC')
hold off

% % figure(8886)
% % subplot(3,1,1)
% % r_max_plot = 6; % Max range to plot [km]
% % [~,ind_km_max] = min(abs(Range.rkm-r_max_plot));
% % imAlpha = ones(size(Alpha.alpha_0));
% % colorLimits = [-1e-4, 3e-4];
% % colorLimits = [0, 3e-4];
% % imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),Alpha.alpha_0(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
% % hold on
% % % plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b')
% % colors = colormap;
% % colors(1,:) = [0 0 0];%set lowest color black
% % colormap(colors);
% % set(gca, 'YDir','normal')
% % set(gca,'color',[1 1 1]);%color background white
% % colorbar
% % %title(sprintf('O2 alpha 0 %s to %s (UTC)',Time.date_ts(1),Time.date_ts(end)))
% % title('Zero order DIAL O_2 absorption')
% % ylabel('Range (km)')
% % title(colorbar,'m^{-1}')
% % timeAxis(Format,Time)
% % %xlim([datenum(datetime(2022,04,30,12,0,0,'TimeZone','UTC')) inf])
% % set(gca, 'XTickLabel', [])
% % hold off
% % subplot(3,1,2)
% % colorLimits = [-5e-5, 5e-5];
% % imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),Alpha.alpha_1(1:ind_km_max,:)+Alpha.alpha_2(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
% % hold on
% % % plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b')
% % colors = colormap;
% % colors(1,:) = [0 0 0];%set lowest color black
% % colormap(colors);
% % set(gca, 'YDir','normal')
% % set(gca,'color',[1 1 1]);%color background white
% % colorbar
% % title(sprintf('Absorption Correction'))
% % ylabel('Range (km)')
% % title(colorbar,'m^{-1}')
% % timeAxis(Format,Time)
% % %xlim([datenum(datetime(2022,04,30,12,0,0,'TimeZone','UTC')) inf])
% % set(gca, 'XTickLabel', [])
% % hold off
% % subplot(3,1,3)
% % colorLimits = [-1e-4, 3e-4];
% % colorLimits = [0, 3e-4];
% % imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),Alpha.alpha_total_raw(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
% % hold on
% % % plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b')
% % colors = colormap;
% % colors(1,:) = [0 0 0];%set lowest color black
% % colormap(colors);
% % set(gca, 'YDir','normal')
% % set(gca,'color',[1 1 1]);%color background white
% % colorbar
% % title('Absorption Final')
% % ylabel('Range (km)')
% % title(colorbar,'m^{-1}')
% % timeAxis(Format,Time)
% % %xlim([datenum(datetime(2022,04,30,12,0,0,'TimeZone','UTC')) inf])
% % hold off



%=WV 2D absorption
% figure(8887)
% subplot(4,1,1)
% r_max_plot = 6; % Max range to plot [km]
% [~,ind_km_max] = min(abs(Range.rkm-r_max_plot));
% imAlpha = ones(size(Alpha.alpha_0wv));
% colorLimits = [-1e-4, 6e-4];
% imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),Alpha.alpha_0wv(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
% hold on
% plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b')
% colors = colormap;
% colors(1,:) = [0 0 0];%set lowest color black
% colormap(colors);
% set(gca, 'YDir','normal')
% set(gca,'color',[1 1 1]);%color background white
% colorbar
% title(sprintf('alpha 0 %s to %s (UTC)',Time.date_ts(1),Time.date_ts(end)))
% ylabel('Range (km)')
% title(colorbar,'m^{-1}')
% timeAxis(Format,Time)
% set(gca, 'XTickLabel', [])
% hold off
% subplot(4,1,2)
% colorLimits = [-8e-6 8e-6];
% imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),Alpha.alpha_1wv(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
% hold on
% plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b')
% colors = colormap;
% colors(1,:) = [0 0 0];%set lowest color black
% colormap(colors);
% set(gca, 'YDir','normal')
% set(gca,'color',[1 1 1]);%color background white
% colorbar
% title(sprintf('alpha 1'))
% ylabel('Range (km)')
% title(colorbar,'m^{-1}')
% timeAxis(Format,Time)
% set(gca, 'XTickLabel', [])
% hold off
% subplot(4,1,3)
% colorLimits = [-8e-7 8e-7];
% imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),Alpha.alpha_2wv(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
% hold on
% plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b')
% colors = colormap;
% colors(1,:) = [0 0 0];%set lowest color black
% colormap(colors);
% set(gca, 'YDir','normal')
% set(gca,'color',[1 1 1]);%color background white
% colorbar
% title(sprintf('alpha 2'))
% title(colorbar,'m^{-1}')
% ylabel('Range (km)')
% timeAxis(Format,Time)
% set(gca, 'XTickLabel', [])
% hold off
% subplot(4,1,4)
% colorLimits = [-1e-4, 6e-4];
% imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),Alpha.alpha_total_rawwv(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits)
% hold on
% plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b')
% colors = colormap;
% colors(1,:) = [0 0 0];%set lowest color black
% colormap(colors);
% set(gca, 'YDir','normal')
% set(gca,'color',[1 1 1]);%color background white
% colorbar
% title('alpha total')
% ylabel('Range (km)')
% title(colorbar,'m^{-1}')
% timeAxis(Format,Time)
% hold off


%=2D BSR
figure(7473)
BSRm = HSRL.BSR;
%imagesc(datenum(Time.date_ts),Range.rkm,BSRm)
h = pcolor(Time.date_ts,Range.rkm,BSRm);
set(h, 'EdgeColor', 'none');
hold on
plot(Time.date_ts([p_point(1) p_point(end)]),[Range.rkm(1) Range.rkm(end)],'--b')
hold off
colorbar
colormap(flipud(hot))
caxis([1 4])
ylim([0 4])
ylabel('Range (m)')
set(gca,'ColorScale','log')
title('BSR')

% figure(74734)
% BSRm = HSRL.BSR;
% imagesc(datenum(Time.date_ts),Range.rkm,HSRL.Ba)
% hold on
% plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b')
% hold off
% colorbar
% colormap(flipud(hot))
% set(gca,'ColorScale','log')
% caxis([10^-8 10^-3])
% timeAxis(Format,Time)
% ylabel('Range (m)')
% set(gca, 'YDir','normal')
% title('BSR')

% figure(7474)
% BSRm = HSRL.BSR;
% imagesc(datenum(Time.date_ts),Range.rkm,HSRL.BSR828-HSRL.BSR)
% hold on
% plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b')
% hold off
% colorbar
% colormap(redblue(64))
% %set(gca,'ColorScale','log')
% caxis([-5 5])
% timeAxis(Format,Time)
% ylabel('Range (m)')
% set(gca, 'YDir','normal')
% title('BSR')


%==2D comparison retrieved vs model
% figure(74744)
% imagesc(datenum(Time.date_ts),Range.rkm,Temperature.T_finalm-Model.T)
% hold on
% plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b')
% hold off
% colorbar
% colormap(redblue(64))
% %set(gca,'ColorScale','log')
% caxis([-20 20])
% timeAxis(Format,Time)
% ylabel('Range (m)')
% set(gca, 'YDir','normal')
% title('T final - T model (K)')


%=T Sonde and DIAL scatter plot
figure(67755)
[~,sondeCut]=min(abs(Range.rm-150));
BL = -20+273;
UR = 25+273;
tempProb = zeros(UR-BL+1);
for ii = 1:size(Sonde.sonde_ind,2)
    if ii==12
    else
    for i = 1:length(BL:UR)
        for j = 1:length(BL:UR)
            tempProb(i,j) = tempProb(i,j)+nnz((diag(Temperature.T_finalm(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))<(BL+i-1)+1)&(diag(Temperature.T_finalm(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))>=(BL+i-1)-0)& (Sonde.T_sonde(sondeCut:end,ii)<(BL+j-1)+1)&(Sonde.T_sonde(sondeCut:end,ii)>=(BL+j-1)-0));
        end
    end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc((BL:UR)-273,(BL:UR)-273,tempProb','AlphaData',imAlpha')
colormap(flipud(hot))%colormap hot
%set(gca,'ColorScale','log')
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
hold on
plot((BL:UR)-273,(BL:UR)-273,'k--')
plot((BL:UR)-273-2,(BL:UR)-273,'k')
plot((BL:UR)-273+2,(BL:UR)-273,'k')
grid on
hold off
set(gca, 'YDir','normal')
xlabel('Temp retrieval (^oC)')
ylabel('Sonde (^oC)')

%=WV Sonde and DIAL scatter plot
% % figure(6775551)
% % [~,sondeCut]=min(abs(Range.rm-150));
% % minEdge =0;
% % maxEdge = 12;
% % binsPer = 4;
% % maxEdge = maxEdge*binsPer;
% % tempProb = zeros(maxEdge-minEdge+1);
% % for ii = 1:size(Sonde.sonde_ind,2)
% %     for i = 1:length(minEdge:maxEdge)
% %         for j = 1:length(minEdge:maxEdge)
% %             tempProb(i,j) = tempProb(i,j)+nnz((diag(AbsHumm(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))<(minEdge+i-1)/binsPer+1)&(diag(AbsHumm(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))>=(minEdge+i-1)/binsPer-0)& (Sonde.AbsHum(sondeCut:end,ii)<(minEdge+j-1)/binsPer+1)&(Sonde.AbsHum(sondeCut:end,ii)>=(minEdge+j-1)/binsPer-0));
% %         end
% %     end
% % end
% % tempProb(tempProb==0)=nan;
% % absHum = ones(Range.i_range,size(Sonde.sonde_ind,2));
% % for iii = 1:size(Sonde.sonde_ind,2)
% %     absHum(:,iii) = diag(AbsHumm(:,Sonde.sonde_ind(:,iii)));
% % end
% % 
% % if ~isempty(Sonde.sonde_ind)
% % absHumReal = absHum(~isnan(absHum) & absHum>0 & Sonde.AbsHum>0);
% % sondeAbsHumReal = Sonde.AbsHum(~isnan(absHum) & absHum>0 & Sonde.AbsHum>0);
% % 
% %     if ~isnan(absHumReal)
% %     FO = fit(absHumReal(:),sondeAbsHumReal(:),'poly1');
% %     else
% %             FO.p1 = nan;
% %         FO.p2 = nan;
% %     end
% % else
% %     FO.p1 = nan;
% %     FO.p2 = nan;
% % end
% % 
% % tempProb = tempProb';
% % tempProb(tempProb==0)=nan;
% % imAlpha=ones(size(tempProb));
% % imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
% % imagesc((minEdge:maxEdge)/binsPer,(minEdge:maxEdge)/binsPer,tempProb','AlphaData',imAlpha')
% % set(gca,'Color','#D3D3D3')
% % a = colorbar;
% % a.Label.String = 'Occurrences';
% % hold on
% % plot((minEdge:maxEdge)/binsPer,(minEdge:maxEdge)/binsPer,'k','LineWidth',2)
% % if ~isnan(FO.p1)
% % plot(FO)
% % end
% % grid on
% % hold off
% % set(gca, 'YDir','normal')
% % title(sprintf('WV total correction\n fit a=%f b=%f',FO.p1,FO.p2))
% % xlabel('WV retrieval (g/m^3)')
% % ylabel('Sonde (g/m^3)')

%=WV Sonde and DIAL scatter plot
% % figure(6775552)
% % [~,sondeCut]=min(abs(Range.rm-150));
% % minEdge =0;
% % maxEdge = 12;
% % binsPer = 4;
% % maxEdge = maxEdge*binsPer;
% % tempProb = zeros(maxEdge-minEdge+1);
% % for ii = 1:size(Sonde.sonde_ind,2)
% %     for i = 1:length(minEdge:maxEdge)
% %         for j = 1:length(minEdge:maxEdge)
% %             tempProb(i,j) = tempProb(i,j)+nnz((diag(AbsHum0m(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))<(minEdge+i-1)/binsPer+1)&(diag(AbsHum0m(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))>=(minEdge+i-1)/binsPer-0)& (Sonde.AbsHum(sondeCut:end,ii)<(minEdge+j-1)/binsPer+1)&(Sonde.AbsHum(sondeCut:end,ii)>=(minEdge+j-1)/binsPer-0));
% %         end
% %     end
% % end
% % tempProb(tempProb==0)=nan;
% % imAlpha=ones(size(tempProb));
% % imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
% % imagesc((minEdge:maxEdge)/binsPer,(minEdge:maxEdge)/binsPer,tempProb,'AlphaData',imAlpha)
% % set(gca,'Color','#D3D3D3')
% % a = colorbar;
% % a.Label.String = 'Occurrences';
% % hold on
% % plot((minEdge:maxEdge)/binsPer,(minEdge:maxEdge)/binsPer,'k','LineWidth',2)
% % grid on
% % hold off
% % set(gca, 'YDir','normal')
% % title('WV 0')
% % xlabel('WV retrieval (g/m^3)')
% % ylabel('Sonde (g/m^3)')

%== scatter plot of model and retrieved
figure(67757)
tempProb = zeros(303-278+1);
for ii = 1:size(Sonde.sonde_ind,2)
    for i = 1:length(278:303)
        for j = 1:length(278:303)
            tempProb(i,j) = tempProb(i,j)+nnz((diag(Model.T(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))<(278+i-1)+1)&(diag(Temperature.T_finalm(sondeCut:end,Sonde.sonde_ind(sondeCut:end,ii)))>=(278+i-1)-0)& (Sonde.T_sonde(sondeCut:end,ii)<(278+j-1)+1)&(Sonde.T_sonde(sondeCut:end,ii)>=(278+j-1)-0));
        end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc((278:303)-273,(278:303)-273,tempProb,'AlphaData',imAlpha)
colormap(flipud(hot))%colormap hot
set(gca,'ColorScale','log')
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
hold on
plot((278:303)-273,(278:303)-273,'k--')
plot((278:303)-273-2,(278:303)-273,'k')
plot((278:303)-273+2,(278:303)-273,'k')
grid on
hold off
set(gca, 'YDir','normal')
xlabel('lapse rate (K)')
ylabel('Sonde (K)')

%=T Sonde-DIAL histogram
figure(67756)
%subplot(1,3,1)
hist11 = nan(1,size(Sonde.sonde_ind,2)*(size(Sonde.sonde_ind,1)-sondeCut));
int = 1;
for ii = 1:size(Sonde.sonde_ind,2)
    if ii==12
    else
    for jj = sondeCut:size(Sonde.sonde_ind,1)
        hist11(int)=Temperature.T_finalm(jj,Sonde.sonde_ind(jj,ii))-Sonde.T_sonde(jj,ii);
        int=int+1;
    end
    end
end
if ~isempty(Sonde.sonde_ind)
    meanHist = mean(hist11,'omitnan');
    stdHist = std(hist11,'omitnan');
    histogram(hist11,'BinWidth',1)
    title(sprintf('Histogram T_{DIAL}-T_{sonde}\n Mean %f, std %f',meanHist,stdHist))
end
xlabel('\DeltaT (^oC)')
xlim([-10 10])
ylabel('Occurrences')

% % subplot(1,3,2)
% % hist11 = nan(1,size(Sonde.sonde_ind,2)*(size(Sonde.sonde_ind,1)-sondeCut));
% % int = 1;
% % for ii = 1:size(Sonde.sonde_ind,2)
% %     if ii==12
% %     elseif Counts.bg_o2off_mol(1,Sonde.sonde_ind(1,ii))>10
% %     for jj = sondeCut:size(Sonde.sonde_ind,1)
% %         hist11(int)=Temperature.T_finalm(jj,Sonde.sonde_ind(jj,ii))-Sonde.T_sonde(jj,ii);
% %         int=int+1;
% %     end
% %     end
% % end
% % if ~isempty(Sonde.sonde_ind)
% %     meanHist = mean(hist11,'omitnan');
% %     stdHist = std(hist11,'omitnan');
% %     histogram(hist11,'BinWidth',1)
% %     title(sprintf('dayHistogram T_{DIAL}-T_{sonde}\n Mean %f, std %f',meanHist,stdHist))
% % end
% % xlabel('\DeltaT (^oC)')
% % xlim([-10 10])
% % ylabel('Occurrences')
% % subplot(1,3,3)
% % hist11 = nan(1,size(Sonde.sonde_ind,2)*(size(Sonde.sonde_ind,1)-sondeCut));
% % int = 1;
% % for ii = 1:size(Sonde.sonde_ind,2)
% %     if ii==12
% %     elseif Counts.bg_o2off_mol(1,Sonde.sonde_ind(1,ii))<=10
% %     for jj = sondeCut:size(Sonde.sonde_ind,1)
% %         hist11(int)=Temperature.T_finalm(jj,Sonde.sonde_ind(jj,ii))-Sonde.T_sonde(jj,ii);
% %         int=int+1;
% %     end
% %     end
% % end
% % if ~isempty(Sonde.sonde_ind)
% %     meanHist = mean(hist11,'omitnan');
% %     stdHist = std(hist11,'omitnan');
% %     histogram(hist11,'BinWidth',1)
% %     title(sprintf('nightHistogram T_{DIAL}-T_{sonde}\n Mean %f, std %f',meanHist,stdHist))
% % end
% % xlabel('\DeltaT (^oC)')
% % xlim([-10 10])
% % ylabel('Occurrences')

%===DIAL and model histogram
figure(67758)
subplot(1,3,1)
hist11 = zeros(1,size(Sonde.sonde_ind,2)*(size(Sonde.sonde_ind,1)-sondeCut));
int = 1;
ModelTm = Model.T;
ModelTm = -6.5/1000.*Range.rm+Model.Ts;
ModelTm(cloud_SDm_above) = NaN;          % Replace mask with NaNs
for ii = 1:size(Sonde.sonde_ind,2)
    for jj = sondeCut:size(Sonde.sonde_ind,1)
        hist11(int)=Sonde.T_sonde(jj,ii)-ModelTm(jj,Sonde.sonde_ind(jj,ii));
        int=int+1;
    end
end
if ~isempty(Sonde.sonde_ind)
    meanHist = mean(hist11,'omitnan');
    stdHist = std(hist11,'omitnan');
    histogram(hist11,'BinWidth',1)
    title(sprintf('Histogram T_{sonde}-T_{Lapse}\n Mean %f, std %f',meanHist,stdHist))
end
xlim([-10 10])
subplot(1,3,2)
%hist11 = zeros(1,size(Sonde.sonde_ind,2)*(size(Sonde.sonde_ind,1)-sondeCut));
hist11 = nan(1,size(Sonde.sonde_ind,2)*(size(Sonde.sonde_ind,1)-sondeCut));
int = 1;
% ModelTm = Model.T;
% ModelTm(cloud_SDm_above) = NaN;          % Replace mask with NaNs
for ii = 1:size(Sonde.sonde_ind,2)
    if Counts.bg_o2off_mol(1,Sonde.sonde_ind(1,ii))>10
    for jj = sondeCut:size(Sonde.sonde_ind,1)
        hist11(int)=Sonde.T_sonde(jj,ii)-ModelTm(jj,Sonde.sonde_ind(jj,ii));
        int=int+1;
    end
    end
end
if ~isempty(Sonde.sonde_ind)
    meanHist = mean(hist11,'omitnan');
    stdHist = std(hist11,'omitnan');
    histogram(hist11,'BinWidth',1)
    title(sprintf('day Histogram T_{sonde}-T_{Lapse}\n Mean %f, std %f',meanHist,stdHist))
end
xlim([-10 10])
subplot(1,3,3)
%hist11 = zeros(1,size(Sonde.sonde_ind,2)*(size(Sonde.sonde_ind,1)-sondeCut));
hist11 = nan(1,size(Sonde.sonde_ind,2)*(size(Sonde.sonde_ind,1)-sondeCut));
int = 1;
% ModelTm = Model.T;
% ModelTm(cloud_SDm_above) = NaN;          % Replace mask with NaNs
for ii = 1:size(Sonde.sonde_ind,2)
    if Counts.bg_o2off_mol(1,Sonde.sonde_ind(1,ii))<=10
    for jj = sondeCut:size(Sonde.sonde_ind,1)
        hist11(int)=Sonde.T_sonde(jj,ii)-ModelTm(jj,Sonde.sonde_ind(jj,ii));
        int=int+1;
    end
    end
end
if ~isempty(Sonde.sonde_ind)
    meanHist = mean(hist11,'omitnan');
    stdHist = std(hist11,'omitnan');
    histogram(hist11,'BinWidth',1)
    title(sprintf('night Histogram T_{sonde}-T_{Lapse}\n Mean %f, std %f',meanHist,stdHist))
end
xlim([-10 10])

%==wv and sonde histogram
% % figure(677564)
% % subplot(2,1,1)
% % hist11 = zeros(1,size(Sonde.sonde_ind,2)*(size(Sonde.sonde_ind,1)-sondeCut));
% % int = 1;
% % for ii = 1:size(Sonde.sonde_ind,2)
% %     for jj = sondeCut:size(Sonde.sonde_ind,1)
% %         hist11(int)=(N_wvm(jj,Sonde.sonde_ind(jj,ii))-Sonde.WV_sonde(jj,ii))./Sonde.WV_sonde(jj,ii)*100;
% %         int=int+1;
% %     end
% % end
% % if ~isempty(Sonde.sonde_ind)
% %     meanHist = mean(hist11,'omitnan');
% %     stdHist = std(hist11,'omitnan');
% %     histogram(hist11,'BinWidth',5)
% %     title(sprintf('Histogram WV_{DIAL}-WV_{sonde} %%\n Mean %f, std %f',meanHist,stdHist))
% % end
% % xlim([-40 40])
% % ylim([0 200])
% % ylabel('Occurrences')
% % 
% % subplot(2,1,2)
% % hist11 = zeros(1,size(Sonde.sonde_ind,2)*(size(Sonde.sonde_ind,1)-sondeCut));
% % int = 1;
% % N_wv0m(N_wv0m==inf | N_wv0m==-inf)=nan;
% % 
% % for ii = 1:size(Sonde.sonde_ind,2)
% %     for jj = sondeCut:size(Sonde.sonde_ind,1)
% %         hist11(int)=(N_wv0m(jj,Sonde.sonde_ind(jj,ii))-Sonde.WV_sonde(jj,ii))./Sonde.WV_sonde(jj,ii)*100;
% %         int=int+1;
% %     end
% % end
% % if ~isempty(Sonde.sonde_ind)
% %     meanHist = mean(hist11,'omitnan');
% %     stdHist = std(hist11,'omitnan');
% %     histogram(hist11,'BinWidth',5)
% %     title(sprintf('Histogram WV0_{DIAL}-WV_{sonde}%%\n Mean %f, std %f',meanHist,stdHist))
% % end
% % %xlabel('\Delta T Sonde-MPD (^oC)')
% % xlim([-40 40])
% % ylim([0 200])
% % ylabel('Occurrences')

%=FFT of Counts
% figure(7726)
% subplot(2,2,1)
% Fs = 1/Range.rangeBin;
% T = Range.rangeBin;
% L = length(Range.rm);
% f = Fs*(0:(L/2))/L;
% Y = fft(Counts.o2on(:,p_point(1)).*Range.rm.^2);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold on
% Y = fft(Counts.o2off(:,p_point(1)).*Range.rm.^2);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% Y = fft(Counts.o2on_mol(:,p_point(1)).*Range.rm.^2);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% Y = fft(Counts.o2off_mol(:,p_point(1)).*Range.rm.^2);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold off
% title('FFT of counts in range at plot point')
% legend('on','off','on mol','off mol')
% xlabel('f (1/m)')
% ylabel('|P1(f)|')
% %figure(722897)
% subplot(2,2,2)
% Y = fft(fillmissing(Alpha.alpha_total(:,p_point(1)),'nearest'));
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold on
% Y = fft(fillmissing(Alpha.alpha_total_raw(:,p_point(1)),'nearest'));
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% Y = fft(fillmissing(Sonde.absorption_sonde{sonde_index}(:,:),'nearest'));
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold off
% legend('Alpha total','alpha total raw','alpha sonde')
% title('FFT of absorption in range dimension at plot point')
% xlabel('Freqency (1/m)')
% %figure(722898)
% subplot(2,2,3)
% Y = fft(fillmissing(Temperature.T_final_tests(:,p_point(1)),'nearest'));
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold on
% Y = fft(fillmissing(Temperature.T_final_test(:,p_point(1)),'nearest'));
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% Y = fft(fillmissing(Sonde.T_sonde(:,sonde_index),'nearest'));
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold off
% legend('T total','T total raw','T sonde')
% title('FFT of Temperature in range dimension at plot point')
% xlabel('Freqency (1/m)')
% %figure(722899)
% subplot(2,2,4)
% Y = fft([ones(1,8) zeros(1,Range.i_range-8)]);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold on
% Y = fft([ones(1,20) zeros(1,Range.i_range-20)]);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% semilogy(f,P1)
% hold off
% legend('8 length','10 length')

%=N bins over time
figure(66572)
bar(Time.date_ts,Counts.NBins)
ylabel('# of bins summed')
%timeAxis(Format,Time)

%=Counts and corrections
figure(6)
plot(Range.rkm,diag(Counts.o2on(:,p_point)),'--','linewidth',2)
hold on
plot(Range.rkm,diag(Counts.o2off(:,p_point)),'--','linewidth',2)
plot(Range.rkm,diag(Counts.o2on_mol(:,p_point)),'linewidth',2)
plot(Range.rkm,diag(Counts.o2off_mol(:,p_point)),'linewidth',2)
%plot(Range.rkm,diag(Counts.o2off_mol(:,p_point)))
%plot(Range.rkm,diag(Counts.wvon(:,p_point)),'--')
%plot(Range.rkm,diag(Counts.wvoff(:,p_point)),'--')
%plot(Range.rm,diag(Counts.delta(:,p_point)),'.-')
%plot(Range.rm,diag(Counts.deltaModel(:,p_point)),'.-')
%plot(Range.rm,diag(Counts.o2on(:,p_point)-Counts.delta(:,p_point)),'--')
%plot(Range.rm,diag(Counts.o2off(:,p_point)-Counts.delta(:,p_point)),'--')
ylim([-200 1000])

ylim([0 4e5])
%load('AfterPulse.mat')
%plot(Range.rm,pulseON(4:160+3))
%plot(Range.rm,[pulseON(4:160); zeros(3,1)])
% plot(Range.rm,diag(Counts.o2on(:,p_point)-pulseON(4:133+3)),'.-')
% plot(Range.rm,diag(Counts.o2off(:,p_point)-pulseON(4:133+3)),'.-')
%%plot(Range.rm,diag(Counts.o2on(:,p_point)-[pulseON(4:160); zeros(3,1)]),'.-')
%%plot(Range.rm,diag(Counts.o2off(:,p_point)-[pulseON(4:160); zeros(3,1)]),'.-')
% plot(Range.rm,diag(Counts.o2on_mol(:,p_point)))
% plot(Range.rm,diag(Counts.o2off_mol(:,p_point)))
%plot(Range.rm(1:end-1),-diff(Counts.o2off(:,p_point(1))))
%%plot(Range.rm(1:end-1),-diff(Counts.o2off(:,p_point(1)))-[pulseON(4:159); zeros(3,1)])

% plot(Range.rm,Model.N_on(:,p_point(1)))
% plot(Range.rm,Model.N_off(:,p_point(1)))
% plot(Range.rm,Model.N_on_pulse(:,p_point(1)),'--')
% plot(Range.rm,Model.N_off_pulse(:,p_point(1)),'--')
grid on
hold off
xlabel('Range (km)')
%legend('On','off','\Delta','On -\Delta','Off -\Delta','Measured afterpulse')
legend('On Com','Off Com','On Mol','Off Col')
title(sprintf(['Counts\n' datestr(Time.date_ts(p_point(1)))]))

% figure(666)
% plot(Range.rm,Model.OverlapOn(:,p_point(1)));
% hold on
% plot(Range.rm,Model.OverlapOff(:,p_point(1)));
% plot(Range.rm,Model.OverlapOn_pulse(:,p_point(1)),'--');
% plot(Range.rm,Model.OverlapOff_pulse(:,p_point(1)),'--');
% hold off
% title('Model counts overlap')
% legend('on','off','on pulse','off pulse')

%=Wavelength overlap correction model
% figure(555)
% plot(Range.rkm,Model.O_on_O_off(:,p_point(1)-10:p_point(1)+10))
% grid on
% title('O_{on}(r)/O_{off}(r) correction')
% xlabel('Range km')

%=Transmission
% % figure(90)
% % subplot(2,1,1)
% % %lnOonOoff = log(Range.rm.^2.*Counts.o2on)-log(Range.rm.^2.*Counts.o2off)-log(Model.transmission.^2);
% % lnOonOoff = log(Range.rm.^2.*Counts.o2on)-log(Range.rm.^2.*Counts.o2off)-log(Sonde.trasmission_sonde{sonde_index}.^2);
% % plot(Range.rm,lnOonOoff(:,p_point(1)))
% % hold on
% % 
% % logicalExc=false(Range.i_range,Time.i_time);
% % alpha_fit_slope = ones(1,Time.i_time);
% % alpha_fit_int=ones(1,Time.i_time);
% % alpha_fit = ones(Range.i_range,Time.i_time);
% % for time = 1:length(Time.ts)
% %         % Set fit exclusion zones
% %         %Temperature.exclusion(:,time) = Range.rm<1500 | Range.rm>2500 | SNRm(:,time)==0 | cloud_SDm_above(:,time) ~= 1;
% %         Temperature.exclusion(:,time) = Range.rm<1500 | Range.rm>2500 |  cloud_SDm_above(:,time);
% %         logicalExc(:,time) = ~logical(Temperature.exclusion(:,time)); 
% %         V = [Range.rm(logicalExc(:,time)) ones(length(Range.rm(logicalExc(:,time))),1)];
% %         %calculate polynomial coefficients     
% %         p = V\Alpha.alpha_0(logicalExc(:,time),time);
% %         alpha_fit_slope(:,time) = p(1);
% %         alpha_fit_int(:,time) = p(2);
% %         alpha_fit(:,time) = alpha_fit_int(:,time) + Range.rm.*alpha_fit_slope(:,time);
% % end
% % transmission = exp(-cumtrapz(Range.rm,alpha_fit));
% % lnOonOoff = log(Range.rm.^2.*Counts.o2on)-log(Range.rm.^2.*Counts.o2off)-log(transmission.^2);
% % %lnOonOoff = lnOonOoff - mean(lnOonOoff(40:67,:),1);
% % lnOonOoff = lnOonOoff - mean(lnOonOoff(14:27,:),1);
% % plot(Range.rm,lnOonOoff(:,p_point(1)),'--')
% % grid on
% % legend('model lnOonOoff','Fit lnOonOoff')
% % hold off
% % 
% % %=log on and off
% % subplot(2,1,2)
% % lnOonOoff=exp(lnOonOoff);
% % plot(Range.rm,lnOonOoff(:,p_point),'--')
% % %plot(Range.rm,lnOonOoff(:,p_point)./Corr(:,p_point))
% % title('lnOonOoff')
% % grid on
% % hold off


%==2D of overlap correction
% figure(9090)
% imagesc(Time.thr,Range.rm/1000,lnOonOoff)
% hold on
% xline(Time.thr(p_point(1)),'K');
% hold off
% set(gca, 'YDir','normal')
% caxis([.9 1.1])
% colorbar
% colormap(redblue(64))
% title('Counts.delta')
% ylabel('Range (km)')
% xlabel('Time (UTC hrs)')
           
% figure
% plot(Range.rm,log(Range.rm.^2 .* Counts.o2on(:,p_point(1)))-log(Range.rm.^2 .* Counts.o2off(:,p_point(1))))
% hold on
% plot(Range.rm,log(Model.transmission(:,p_point(1)).^2))
% 
% %TT = Stemp(minJJ) + lapse(minII)/1000 .* Range.rm;                           %[K] (1 x r) Temperature model as a function of r 
% %PP = Spress(minKK) .* (Stemp(minJJ)./TT).^(-5.2199);                       %[atm] (1 x r) Pressure model as a function of r  
% %alpha = absorption_O2_770_model(TT,PP,Spectrum.nu_online(1),Model.WV(:,p_point)); %[m-1] Funcrtion to calculate theoretical absorption
% %Transmission = exp(-cumtrapz(Range.rm,alpha));
% plot(Range.rm,log(Transmission(:,p_point(1)).^2)+const(p_point(1)))
% hold off
% legend('Counts','model','fit')
% 
% figure
% plot(Range.rm,exp(log(Range.rm.^2 .* Counts.o2on(:,p_point(1)))-log(Range.rm.^2 .* Counts.o2off(:,p_point(1)))-log(Transmission(:,p_point(1)).^2)-const(p_point(1))))
% hold on
% plot(Range.rm,exp(log(Range.rm.^2 .* Counts.o2on(:,p_point(1)))-log(Range.rm.^2 .* Counts.o2off(:,p_point(1)))-log(Model.transmission(:,p_point(1)).^2)))
% %plot(Range.rm,smooth(exp(log(Range.rm.^2 .* Counts.o2on(:,p_point(1)))-log(Range.rm.^2 .* Counts.o2off(:,p_point(1)))-log(Transmission(:,p_point(1)).^2)-const(p_point(1))),5))
% plot(Range.rm,overlapCorr(:,p_point))
% hold off
% legend('Fit','model')

%== online offline background difference
if ishandle(64888)
    figure(64888)
    clf(64888)
end
figure(64888)
line(Time.thr,smooth(-Counts.bg_o2on+Counts.bg_o2off,10),'Color','r')
ylim([-5 5])
ax1=gca;
ax1.XColor = 'r';
ax1.YColor = 'r';
xlabel('time (hrs)')
ylabel('-Online Background + offline background')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(Time.thr,smooth(Alpha.alpha_0(4,:),10),'Parent',ax2,'Color','k')
line(Time.thr,smooth(Alpha.alpha_0(9,:),10),'Parent',ax2,'Color','b')
ylim([0 5e-4])
grid on
ylabel('alpha 0')

%===TSOA Power
% % % figure(3927)
% % % plot(Data.Power.O2Online.TimeStamp,Data.Power.O2Online.LaserPower)
% % % hold on
% % % plot(Data.Power.O2Offline.TimeStamp,Data.Power.O2Offline.LaserPower)
% % % xline(Time.ts(p_point(1))/60/60)
% % % hold off
% % % xlim([Time.ts(1) Time.ts(end)]/60/60)
% % % grid on
% % % title('Power Measuremet')
% % % legend('Online','Offline')

%== Potential temperature==
% figure(337)
% subplot(1,2,1)
% Model.T_pot = Model.T.*(0.9869232667160128./Model.P).^(0.286);
% T_pot = Temperature.T_finalm.*(0.9869232667160128./Model.P).^(0.286);
% T_pot = real(T_pot);
% if ~isempty(Sonde.T_sonde)
%     Sonde.T_pot = Sonde.T_sonde.*(0.9869232667160128./Sonde.P_sonde).^(0.286);
% end
% plot(diag(Model.T_pot(:,p_point)),Range.rkm)
% hold on
% plot(Sonde.T_pot(:,sonde_index),Range.rkm)
% plot(diag(T_pot(:,p_point)),Range.rkm)
% hold off
% legend('model','measurement','sonde')
% grid on
% xlabel('Potential Temperature (K)')
% ylabel('Range (Range.rkm)')
% subplot(1,2,2)
% plot(diff(diag(Model.T_pot(:,p_point)))/Range.rangeBin ,Range.rkm(1:end-1))
% hold on
% plot(diff(diag(T_pot(:,p_point)))/Range.rangeBin,Range.rkm(1:end-1))
% plot(diff(Sonde.T_pot(:,sonde_index))/Range.rangeBin,Range.rkm(1:end-1))
% xlim([-.15 .15])
% hold off
% legend('model','measurement','sonde')
% grid on
% xlabel('Potential Temperature lapse rate (K/m)')
% ylabel('Range (Range.rkm)')
% 
% 
% figure(338)
% r_max_plot = 5; % Max range to plot [km]
% [~,ind_km_max] = min(abs(Range.rkm-r_max_plot)); % Max range Index
% r_min_plot = .5; % Max range to plot [km]
% [~,ind_km_min] = min(abs(Range.rkm-r_min_plot)); % Max range Index
% imAlpha=ones(size(Temperature.T_finalm(1:ind_km_max,:))); % initalize transpernecy matrix
% imAlpha(isnan(Temperature.T_finalm(1:ind_km_max,:)))=0; % Set AlphaData to not plot NaNs
% %imAlpha(cloud_SDm(1:ind_km_max,:) ~= 1)=1; %Set transperency so that cloud mask can be seen
% 
% imAlpha(~cloud_SDm(1:ind_km_max,:))=1; %Set transperency so that cloud mask can be seen
% %imAlpha(cloud_SDm_above(1:ind_km_max,:) ~= 1)=1; %Set transperency so that cloud mask can be seen
% colorLimits = [min(Temperature.T_finalm(ind_km_min:ind_km_max,4:end-4)-1,[],'all'), max(Temperature.T_finalm(ind_km_min:ind_km_max,4:end-4),[],'all')]; % Set color limits to only include data -1
% colorLimits =[290 340];
% imagesc(datenum(Time.date_ts),Range.rkm(1:ind_km_max),T_pot(1:ind_km_max,:),'AlphaData',imAlpha,colorLimits) %Plot
% hold on
% %sonde lines
% plot(datenum(Time.date_ts([p_point(1) p_point(end)])),[Range.rkm(1) Range.rkm(end)],'--b','LineWidth',2)
% %yline(.6)
% for ii = 1:size(Sonde.sonde_ind,2)
%     plot(datenum(Time.date_ts([Sonde.sonde_ind(1,ii) Sonde.sonde_ind(end,ii)])),[Range.rkm(1) Range.rkm(end)],'--k')
% end
% colors = colormap; % Set colors to colormap
% colors(1,:) = [0 0 0]; % Set lowest color black to make clouds appear black
% colormap(colors); % Make colormap new colors
% set(gca, 'YDir','normal') % Set y axis increasing
% set(gca,'color',[1 1 1]);% Color background white
% colorbar % Add colorbar
% timeAxis(Format,Time)
% title(sprintf('Potential Temperature retrieval')) %Add title
% xlabel('Time (UTC hours)')
% ylabel('Range (km)')
% title(colorbar,'Temperature (K)') % Add title to colorbar
% hold off

% figure(7824)
% tempComparison = ones(Range.i_range,size(Sonde.sonde_ind,2));
% for jj = 1:size(Sonde.sonde_ind,2)
%     tempComparison(:,jj) = diag(Temperature.T_finalm(:,Sonde.sonde_ind(:,jj)))-Sonde.T_sonde(:,jj);
% end
% if ~isempty(Sonde.sonde_ind)
%  meanTemp = mean(tempComparison,2,'omitnan');
%  stdTemp = std(tempComparison,0,2,'omitnan');
% else
%     tempComparison = nan;
%     meanTemp = nan;
%     stdTemp = nan;
% end
% numberNonNan = sum(~isnan(tempComparison),2);
% %numberNonNan = 1;
% stderTemp = stdTemp./numberNonNan;
% plot(meanTemp,Range.rkm,'r','linewidth',2)
% hold on
% plot(meanTemp+stdTemp,Range.rkm,'--k')
% plot(meanTemp-stdTemp,Range.rkm,'--k')
% plot(meanTemp+stderTemp,Range.rkm,'--b')
% plot(meanTemp-stderTemp,Range.rkm,'--b')
% ylabel('Range (km)')
% xlabel('\DeltaT MPD-Sonde (^oC)')
% legend('Mean','Sandard deviation')
% hold off
% xlim([-20 20])
% ylim([0 4])
% grid on


%==sonde and DIAL differnces
figure(883469)
%subplot(1,3,1)
tempComparison = nan(Range.i_range,size(Sonde.sonde_ind,2));
for jj = 1:size(Sonde.sonde_ind,2)
    if jj==12
    else
    tempComparison(:,jj) = diag(Temperature.T_finalm(:,Sonde.sonde_ind(:,jj)))-Sonde.T_sonde(:,jj);
    end
end
if ~isempty(Sonde.sonde_ind)
 meanTemp = mean(tempComparison,2,'omitnan');
 stdTemp = std(tempComparison,0,2,'omitnan');
else
    meanTemp = nan;
    stdTemp = nan;
end
tempProb = zeros(Range.i_range,121);
for ii = 1:size(Sonde.sonde_ind,2)
    if ii==12
    else
    for i = 1:Range.i_range
        for j = -60:60
            tempDiff = Temperature.T_finalm(i,Sonde.sonde_ind(i,ii))-Sonde.T_sonde(i,ii);
            tempProb(i,j+61) = tempProb(i,j+61) + double((tempDiff<j+1)&&(tempDiff>=j));
        end
    end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc(-60:60,Range.rkm,tempProb,'AlphaData',imAlpha)
hold on
plot(meanTemp,Range.rkm,'b','linewidth',2)
plot(meanTemp+stdTemp,Range.rkm,'--k')
plot(meanTemp-stdTemp,Range.rkm,'--k')
colormap(flipud(hot))%colormap hot
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
xline(0)
xline(-2)
xline(2)
hold off
ylim([0 4])
xlim([-20 20])
set(gca, 'YDir','normal')
xlabel('\DeltaT MPD-Sonde (^oC)')
ylabel('Range (km)')
grid on

% subplot(1,3,2)
% tempComparison = nan(Range.i_range,size(Sonde.sonde_ind,2));
% for jj = 1:size(Sonde.sonde_ind,2)
%     if jj==12
%     elseif Counts.bg_o2off_mol(1,Sonde.sonde_ind(1,jj))>10
%     tempComparison(:,jj) = diag(Temperature.T_finalm(:,Sonde.sonde_ind(:,jj)))-Sonde.T_sonde(:,jj);
%     end
% end
% if ~isempty(Sonde.sonde_ind)
%  meanTemp = mean(tempComparison,2,'omitnan');
%  stdTemp = std(tempComparison,0,2,'omitnan');
% else
%     meanTemp = nan;
%     stdTemp = nan;
% end
% tempProb = zeros(Range.i_range,121);
% for ii = 1:size(Sonde.sonde_ind,2)
%     if ii==12
%     elseif Counts.bg_o2off_mol(1,Sonde.sonde_ind(1,ii))>10
%     for i = 1:Range.i_range
%         for j = -60:60
%             tempDiff = Temperature.T_finalm(i,Sonde.sonde_ind(i,ii))-Sonde.T_sonde(i,ii);
%             tempProb(i,j+61) = tempProb(i,j+61) + double((tempDiff<j+1)&&(tempDiff>=j));
%         end
%     end
%     end
% end
% tempProb(tempProb==0)=nan;
% imAlpha=ones(size(tempProb));
% imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
% imagesc(-60:60,Range.rkm,tempProb,'AlphaData',imAlpha)
% hold on
% plot(meanTemp,Range.rkm,'b','linewidth',2)
% plot(meanTemp+stdTemp,Range.rkm,'--k')
% plot(meanTemp-stdTemp,Range.rkm,'--k')
% colormap(flipud(hot))%colormap hot
% set(gca,'Color','#D3D3D3')
% a = colorbar;
% a.Label.String = 'Occurrences';
% xline(0)
% xline(-2)
% xline(2)
% hold off
% ylim([0 4])
% xlim([-20 20])
% set(gca, 'YDir','normal')
% xlabel('\DeltaT MPD-Sonde (^oC)')
% ylabel('Range (km)')
% grid on
% 
% subplot(1,3,3)
% tempComparison = nan(Range.i_range,size(Sonde.sonde_ind,2));
% for jj = 1:size(Sonde.sonde_ind,2)
%     if jj==12
%     elseif Counts.bg_o2off_mol(1,Sonde.sonde_ind(1,jj))<=10
%     tempComparison(:,jj) = diag(Temperature.T_finalm(:,Sonde.sonde_ind(:,jj)))-Sonde.T_sonde(:,jj);
%     end
% end
% if ~isempty(Sonde.sonde_ind)
%  meanTemp = mean(tempComparison,2,'omitnan');
%  stdTemp = std(tempComparison,0,2,'omitnan');
% else
%     meanTemp = nan;
%     stdTemp = nan;
% end
% tempProb = zeros(Range.i_range,121);
% for ii = 1:size(Sonde.sonde_ind,2)
%     if ii==12
%     elseif Counts.bg_o2off_mol(1,Sonde.sonde_ind(1,ii))<=10
%     for i = 1:Range.i_range
%         for j = -60:60
%             tempDiff = Temperature.T_finalm(i,Sonde.sonde_ind(i,ii))-Sonde.T_sonde(i,ii);
%             tempProb(i,j+61) = tempProb(i,j+61) + double((tempDiff<j+1)&&(tempDiff>=j));
%         end
%     end
%     end
% end
% tempProb(tempProb==0)=nan;
% imAlpha=ones(size(tempProb));
% imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
% imagesc(-60:60,Range.rkm,tempProb,'AlphaData',imAlpha)
% hold on
% plot(meanTemp,Range.rkm,'b','linewidth',2)
% plot(meanTemp+stdTemp,Range.rkm,'--k')
% plot(meanTemp-stdTemp,Range.rkm,'--k')
% colormap(flipud(hot))%colormap hot
% set(gca,'Color','#D3D3D3')
% a = colorbar;
% a.Label.String = 'Occurrences';
% xline(0)
% xline(-2)
% xline(2)
% hold off
% ylim([0 4])
% xlim([-20 20])
% set(gca, 'YDir','normal')
% xlabel('\DeltaT MPD-Sonde (^oC)')
% ylabel('Range (km)')
% grid on



figure(88346999)
%subplot(1,3,1)
tempComparison = nan(Range.i_range,size(Sonde.sonde_ind,2));
for jj = 1:size(Sonde.sonde_ind,2)
    if jj==12
    else
    tempComparison(:,jj) = (diag(Temperature.T_finalm(:,Sonde.sonde_ind(:,jj)))-Sonde.T_sonde(:,jj))./Sonde.T_sonde(:,jj)*100;
    end
end
if ~isempty(Sonde.sonde_ind)
 meanTemp = mean(tempComparison,2,'omitnan');
 stdTemp = std(tempComparison,0,2,'omitnan');
else
    meanTemp = nan;
    stdTemp = nan;
end
tempProb = zeros(Range.i_range,121);
for ii = 1:size(Sonde.sonde_ind,2)
    if ii==12
    else
    for i = 1:Range.i_range
        for j = -60:60
            tempDiff = (Temperature.T_finalm(i,Sonde.sonde_ind(i,ii))-Sonde.T_sonde(i,ii))./Sonde.T_sonde(i,ii)*100;
            tempProb(i,j+61) = tempProb(i,j+61) + double((tempDiff<j+1)&&(tempDiff>=j));
        end
    end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc(-60:60,Range.rkm,tempProb,'AlphaData',imAlpha)
hold on
plot(meanTemp,Range.rkm,'b','linewidth',2)
plot(meanTemp+stdTemp,Range.rkm,'--k')
plot(meanTemp-stdTemp,Range.rkm,'--k')
colormap(flipud(hot))%colormap hot
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
xline(0)
xline(-2)
xline(2)
hold off
ylim([0 4])
xlim([-20 20])
set(gca, 'YDir','normal')
xlabel('\DeltaT MPD-Sonde (%)')
ylabel('Range (km)')
grid on

figure(2345)
plot(Temperature.TempStd(:,p_point(1)),Range.rkm)
hold on
plot(Temperature.TempStds(:,p_point(1)),Range.rkm)
plot(stdTemp,Range.rkm)
hold off
xlabel('errest')
    ylabel('Range (km)')
    legend('Temperature std', 'Temp std smooth','MDP-Sonde std')
    grid on
    


figure(8883469)
%factor = 1e6;
factor = 1;
tempComparison = nan(Range.i_range,size(Sonde.sonde_ind,2));
for jj = 1:size(Sonde.sonde_ind,2)
    if jj==12
    else
    tempComparison(:,jj) = (diag(Alpha.alpha_totalm(:,Sonde.sonde_ind(:,jj)))-Sonde.absorption_sonde{jj})./Sonde.absorption_sonde{jj}*100;
    end
end
if ~isempty(Sonde.sonde_ind)
 meanTemp = mean(tempComparison,2,'omitnan');
 stdTemp = std(tempComparison,0,2,'omitnan');
else
    meanTemp = nan;
    stdTemp = nan;
end
tempProb = zeros(Range.i_range,121);
for ii = 1:size(Sonde.sonde_ind,2)
    if ii==12
    else
    for i = 1:Range.i_range
        for j = -60:60
            tempDiff = (Alpha.alpha_totalm(i,Sonde.sonde_ind(i,ii))-Sonde.absorption_sonde{ii}(i))./Sonde.absorption_sonde{ii}(i)*100;
            tempProb(i,j+61) = tempProb(i,j+61) + double((tempDiff<j+1)&&(tempDiff>=j))/factor;
        end
    end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc(-60:60,Range.rkm,tempProb,'AlphaData',imAlpha)
hold on
plot((meanTemp)*factor,Range.rkm,'b','linewidth',2)
plot((meanTemp+stdTemp)*factor,Range.rkm,'--k')
plot((meanTemp-stdTemp)*factor,Range.rkm,'--k')
colormap(flipud(hot))%colormap hot
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
xline(0)
xline(-2)
xline(2)
hold off
ylim([0 4])
%xlim([-20 20])
set(gca, 'YDir','normal')
xlabel('\Delta a MPD-Sonde (%)')
ylabel('Range (km)')
grid on

% figure(888346999)
% factor = 1e6;
% factor = 1;
% tempComparison = nan(Range.i_range,size(Sonde.sonde_ind,2));
% for jj = 1:size(Sonde.sonde_ind,2)
%     if jj==12
%     else
%     tempComparison(:,jj) = (diag(Alpha.alphatotals(:,Sonde.sonde_ind(:,jj)))-Sonde.absorption_sonde{jj})./Sonde.absorption_sonde{jj}*100;
%     %tempComparison(:,jj) = (diag(Alpha.alphatotals(:,Sonde.sonde_ind(:,jj)))-Sonde.absorption_sonde{jj});
%     end
% end
% if ~isempty(Sonde.sonde_ind)
%  meanTemp = mean(tempComparison,2,'omitnan');
%  stdTemp = std(tempComparison,0,2,'omitnan');
% else
%     meanTemp = nan;
%     stdTemp = nan;
% end
% tempProb = zeros(Range.i_range,121);
% for ii = 1:size(Sonde.sonde_ind,2)
%     if ii==12
%     else
%     for i = 1:Range.i_range
%         for j = -60:60
%             tempDiff = (Alpha.alpha_totalm(i,Sonde.sonde_ind(i,ii))-Sonde.absorption_sonde{ii}(i))./Sonde.absorption_sonde{ii}(i)*100;
%            % tempDiff = (Alpha.alpha_totalm(i,Sonde.sonde_ind(i,ii))-Sonde.absorption_sonde{ii}(i));
%             tempProb(i,j+61) = tempProb(i,j+61) + double((tempDiff<j+1)&&(tempDiff>=j))/factor;
%         end
%     end
%     end
% end
% tempProb(tempProb==0)=nan;
% imAlpha=ones(size(tempProb));
% imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
% imagesc(-60:60,Range.rkm,tempProb,'AlphaData',imAlpha)
% hold on
% plot((meanTemp)*factor,Range.rkm,'b','linewidth',2)
% plot((meanTemp+stdTemp)*factor,Range.rkm,'--k')
% plot((meanTemp-stdTemp)*factor,Range.rkm,'--k')
% colormap(flipud(hot))%colormap hot
% set(gca,'Color','#D3D3D3')
% a = colorbar;
% a.Label.String = 'Occurrences';
% xline(0)
% xline(-2)
% xline(2)
% hold off
% ylim([0 4])
% %xlim([-20 20])
% set(gca, 'YDir','normal')
% xlabel('\Delta a MPD-Sonde (%)')
% ylabel('Range (km)')
% grid on

figure(777)
tempComparison = nan(Range.i_range,size(Sonde.sonde_ind,2));
for jj = 1:size(Sonde.sonde_ind,2)
    if jj==12
    else
    tempComparison(:,jj) = diag(Alpha.alpha_totalm(:,Sonde.sonde_ind(:,jj)))-Sonde.absorption_sonde{jj};
    end
end
if ~isempty(Sonde.sonde_ind)
 meanTemp = mean(tempComparison,2,'omitnan');
 stdTemp = std(tempComparison,0,2,'omitnan');
else
    meanTemp = nan;
    stdTemp = nan;
end
plot(stdTemp,Range.rkm)
hold on
plot(Alpha.alpha_0_err(:,p_point(1)),Range.rkm)
plot(Alpha.alpha_total_err(:,p_point(1)),Range.rkm)
plot(Alpha.alpha_total_errs(:,p_point(1)),Range.rkm)
plot(std(Alpha.alpha_total_rawf(:,p_point(1),:),0,3),Range.rkm)
hold off
legend('Sonde-alpha std','alpha 0 err','total err','total alpha err smooth')
xlabel('absorption diff (1/m)')
ylabel('Range km')


%==sonde dial and model differnces
figure(8834699)
tempComparison = ones(Range.i_range,size(Sonde.sonde_ind,2));
for jj = 1:size(Sonde.sonde_ind,2)
    tempComparison(:,jj) = diag(Temperature.T_finalm(:,Sonde.sonde_ind(:,jj)))-diag(Model.T(:,Sonde.sonde_ind(:,jj)));
end
if ~isempty(Sonde.sonde_ind)
 meanTemp = mean(tempComparison,2,'omitnan');
 stdTemp = std(tempComparison,0,2,'omitnan');
else
    meanTemp = nan;
    stdTemp = nan;
end

subplot(1,2,2)
tempProb = zeros(Range.i_range,121);
for ii = 1:size(Sonde.sonde_ind,2)
    for i = 1:Range.i_range
        for j = -60:60
            tempDiff = Temperature.T_finalm(i,Sonde.sonde_ind(i,ii))-Model.T(i,Sonde.sonde_ind(i,ii));
            tempProb(i,j+61) = tempProb(i,j+61) + double((tempDiff<j+1)&&(tempDiff>=j));
        end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc(-60:60,Range.rkm,tempProb,'AlphaData',imAlpha)
hold on
plot(meanTemp,Range.rkm,'b','linewidth',2)
plot(meanTemp+stdTemp,Range.rkm,'--k')
plot(meanTemp-stdTemp,Range.rkm,'--k')
colormap(flipud(hot))%colormap hot
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
xline(0)
xline(-2)
xline(2)
hold off
ylim([0 4])
xlim([-20 20])
set(gca, 'YDir','normal')
xlabel('\DeltaT MPD-Model (^oC)')
ylabel('Range (km)')
grid on
subplot(1,2,1)

tempComparison = ones(Range.i_range,size(Sonde.sonde_ind,2));
for jj = 1:size(Sonde.sonde_ind,2)
    tempComparison(:,jj) = Sonde.T_sonde(:,jj)-diag(Model.T(:,Sonde.sonde_ind(:,jj)));
end
if ~isempty(Sonde.sonde_ind)
 meanTemp = mean(tempComparison,2,'omitnan');
 stdTemp = std(tempComparison,0,2,'omitnan');
else
    meanTemp = nan;
    stdTemp = nan;
end

tempProb = zeros(Range.i_range,121);
for ii = 1:size(Sonde.sonde_ind,2)
    for i = 1:Range.i_range
        for j = -60:60
            tempDiff = Sonde.T_sonde(i,ii)-Model.T(i,Sonde.sonde_ind(i,ii));
            tempProb(i,j+61) = tempProb(i,j+61) + double((tempDiff<j+1)&&(tempDiff>=j));
        end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc(-60:60,Range.rkm,tempProb,'AlphaData',imAlpha)
hold on
plot(meanTemp,Range.rkm,'b','linewidth',2)
plot(meanTemp+stdTemp,Range.rkm,'--k')
plot(meanTemp-stdTemp,Range.rkm,'--k')
colormap(flipud(hot))%colormap hot
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
xline(0)
xline(-2)
xline(2)
hold off
ylim([0 4])
xlim([-20 20])
set(gca, 'YDir','normal')
xlabel('\DeltaT Sonde-Model (^oC)')
ylabel('Range (km)')
grid on


figure(11111)
subplot(1,2,1)

fitTemp =(Temperature.L_fit_sm_test(:,:,end).*Range.rm+Temperature.Ts_fit(:,:,end));
tempComparison = ones(Range.i_range,size(Sonde.sonde_ind,2));
for jj = 1:size(Sonde.sonde_ind,2)
    tempComparison(:,jj) = Sonde.T_sonde(:,jj)-diag(fitTemp(:,Sonde.sonde_ind(:,jj)));
end
if ~isempty(Sonde.sonde_ind)
 meanTemp = mean(tempComparison,2,'omitnan');
 stdTemp = std(tempComparison,0,2,'omitnan');
else
    meanTemp = nan;
    stdTemp = nan;
end

tempProb = zeros(Range.i_range,121);
for ii = 1:size(Sonde.sonde_ind,2)
    for i = 1:Range.i_range
        for j = -60:60
            tempDiff = Sonde.T_sonde(i,ii)-fitTemp(i,Sonde.sonde_ind(i,ii));
            tempProb(i,j+61) = tempProb(i,j+61) + double((tempDiff<j+1)&&(tempDiff>=j));
        end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc(-60:60,Range.rkm,tempProb,'AlphaData',imAlpha)
hold on
plot(meanTemp,Range.rkm,'b','linewidth',2)
plot(meanTemp+stdTemp,Range.rkm,'--k')
plot(meanTemp-stdTemp,Range.rkm,'--k')
colormap(flipud(hot))%colormap hot
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
xline(0)
xline(-2)
xline(2)
hold off
ylim([0 4])
xlim([-20 20])
set(gca, 'YDir','normal')
xlabel('\DeltaT Sonde-fit (^oC)')
ylabel('Range (km)')
grid on

subplot(1,2,2)

tempComparison = ones(Range.i_range,size(Sonde.sonde_ind,2));
for jj = 1:size(Sonde.sonde_ind,2)
    tempComparison(:,jj) = diag(Temperature.T_finalm(:,Sonde.sonde_ind(:,jj)))-diag(fitTemp(:,Sonde.sonde_ind(:,jj)));
end
if ~isempty(Sonde.sonde_ind)
 meanTemp = mean(tempComparison,2,'omitnan');
 stdTemp = std(tempComparison,0,2,'omitnan');
else
    meanTemp = nan;
    stdTemp = nan;
end

tempProb = zeros(Range.i_range,121);
for ii = 1:size(Sonde.sonde_ind,2)
    for i = 1:Range.i_range
        for j = -60:60
            tempDiff = Temperature.T_finalm(i,Sonde.sonde_ind(i,ii))-fitTemp(i,Sonde.sonde_ind(i,ii));
            tempProb(i,j+61) = tempProb(i,j+61) + double((tempDiff<j+1)&&(tempDiff>=j));
        end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc(-60:60,Range.rkm,tempProb,'AlphaData',imAlpha)
hold on
plot(meanTemp,Range.rkm,'b','linewidth',2)
plot(meanTemp+stdTemp,Range.rkm,'--k')
plot(meanTemp-stdTemp,Range.rkm,'--k')
colormap(flipud(hot))%colormap hot
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
xline(0)
xline(-2)
xline(2)
hold off
ylim([0 4])
xlim([-20 20])
set(gca, 'YDir','normal')
xlabel('\DeltaT MPD-fit (^oC)')
ylabel('Range (km)')
grid on

figure(1112)

tempComparison = Temperature.T_finalm-fitTemp;
 meanTemp = mean(tempComparison,2,'omitnan');
 stdTemp = std(tempComparison,0,2,'omitnan');

 tempProb = zeros(Range.i_range,121);
for ii = 1:Time.i_time
    for i = 1:Range.i_range
        for j = -60:60
            tempDiff = Temperature.T_finalm(i,ii)-fitTemp(i,ii);
            tempProb(i,j+61) = tempProb(i,j+61) + double((tempDiff<j+1)&&(tempDiff>=j));
        end
    end
end
tempProb(tempProb==0)=nan;
imAlpha=ones(size(tempProb));
imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
imagesc(-60:60,Range.rkm,tempProb,'AlphaData',imAlpha)
hold on
plot(meanTemp,Range.rkm,'b','linewidth',2)
plot(meanTemp+stdTemp,Range.rkm,'--k')
plot(meanTemp-stdTemp,Range.rkm,'--k')
colormap(flipud(hot))%colormap hot
set(gca,'Color','#D3D3D3')
a = colorbar;
a.Label.String = 'Occurrences';
xline(0)
xline(-2)
xline(2)
hold off
ylim([0 4])
xlim([-20 20])
set(gca, 'YDir','normal')
xlabel('\DeltaT Temperature.T_finalm-fit (^oC)')
ylabel('Range (km)')
grid on


%==Sonde and DIAL absorption differnces
% % figure(78245)
% % tempComparison = ones(Range.i_range,size(Sonde.sonde_ind,2));
% % for jj = 1:size(Sonde.sonde_ind,2)
% %     tempComparison(:,jj) = (diag(Alpha.alpha_totalm(:,Sonde.sonde_ind(:,jj)))-Sonde.absorption_sonde{jj})./Sonde.absorption_sonde{jj};
% %     sondeabsorption(:,jj)=Sonde.absorption_sonde{jj};
% % end
% % if ~isempty(Sonde.sonde_ind)
% % meanTemp = mean(tempComparison,2,'omitnan');
% % stdTemp = std(tempComparison,0,2,'omitnan');
% % else
% %     tempComparison = nan;
% %     meanTemp = nan;
% %     stdTemp = nan;
% % end
% % numberNonNan = sum(~isnan(tempComparison),2);
% % %numberNonNan = 1;
% % stderTemp = stdTemp./numberNonNan;
% % 
% % for jj = 1:size(Sonde.sonde_ind,2)
% % plot((diag(Alpha.alpha_totalm(:,Sonde.sonde_ind(:,jj)))-sondeabsorption(:,jj))./sondeabsorption(:,jj),Range.rkm)
% % hold on
% % end
% % 
% % plot(meanTemp,Range.rkm,'r','linewidth',2)
% % plot(meanTemp+stdTemp,Range.rkm,'--k')
% % plot(meanTemp-stdTemp,Range.rkm,'--k')
% % plot(meanTemp+stderTemp,Range.rkm,'--b')
% % plot(meanTemp-stderTemp,Range.rkm,'--b')
% % 
% % 
% % ylabel('Range (km)')
% % xlabel('\Delta absorption %')
% % legend('Mean','Sandard deviation')
% % hold off
% % xlim([-3 3])
% % ylim([0 4])
% % grid on

%===wv and sonde differences
% % figure(883470)
% % tempComparison = ones(Range.i_range,size(Sonde.sonde_ind,2));
% % divider = 1;
% % tempProb = zeros(Range.i_range,121);
% % for ii = 1:size(Sonde.sonde_ind,2)
% %     tempComparison(:,ii) = (diag(N_wv0m(:,Sonde.sonde_ind(:,ii))-Sonde.WV_sonde(:,ii))./Sonde.WV_sonde(:,ii))*100;
% %     for i = 1:Range.i_range
% %         for j = -60:60
% %             tempDiff = ((N_wv0m(i,Sonde.sonde_ind(i,ii))-Sonde.WV_sonde(i,ii))./Sonde.WV_sonde(i,ii))*100;
% %             tempProb(i,j+61) = tempProb(i,j+61) + double((tempDiff<(j+1)/divider)&&(tempDiff>=j/divider));
% %         end
% %     end
% % end
% % if ~isempty(Sonde.sonde_ind)
% % meanTemp = mean(tempComparison,2,'omitnan');
% % stdTemp = std(tempComparison,0,2,'omitnan');
% % else
% %     tempComparison = nan;
% %     meanTemp = nan;
% %     stdTemp = nan;
% % end
% % tempProb(tempProb==0)=nan;
% % imAlpha=ones(size(tempProb));
% % imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
% % imagesc(-60:60/divider,Range.rkm,tempProb,'AlphaData',imAlpha)
% % hold on
% % plot(meanTemp,Range.rkm,'b','linewidth',2)
% % plot(meanTemp+stdTemp,Range.rkm,'--k')
% % plot(meanTemp-stdTemp,Range.rkm,'--k')
% % hold off
% % colormap(flipud(hot))%colormap hot
% % set(gca,'Color','#D3D3D3')
% % a = colorbar;
% % a.Label.String = 'Occurrences';
% % hold on
% % xline(0)
% % hold off
% % ylim([0 4])
% % xlim([-30 30])
% % set(gca, 'YDir','normal')
% % xlabel('\Delta WV % MPD-Sonde (^oC)')
% % ylabel('Range (km)')
% % title('WV0')
% % grid on


%==sonde and WV differnces number density
% figure(883471)
% tempComparison = ones(Range.i_range,size(Sonde.sonde_ind,2));
% divider = 1;
% tempProb = zeros(Range.i_range,121);
% for ii = 1:size(Sonde.sonde_ind,2)
%     tempComparison(:,ii) = (diag(N_wvm(:,Sonde.sonde_ind(:,ii))-Sonde.WV_sonde(:,ii))./Sonde.WV_sonde(:,ii))*100;
%     for i = 1:Range.i_range
%         for j = -60:60
%             tempDiff = ((N_wvm(i,Sonde.sonde_ind(i,ii))-Sonde.WV_sonde(i,ii))./Sonde.WV_sonde(i,ii))*100;
%             tempProb(i,j+61) = tempProb(i,j+61) + double((tempDiff<(j+1)/divider)&&(tempDiff>=j/divider));
%         end
%     end
% end
% if ~isempty(Sonde.sonde_ind)
% meanTemp = mean(tempComparison,2,'omitnan');
% stdTemp = std(tempComparison,0,2,'omitnan');
% else
%     tempComparison = nan;
%     meanTemp = nan;
%     stdTemp = nan;
% end
% tempProb(tempProb==0)=nan;
% imAlpha=ones(size(tempProb));
% imAlpha(isnan(tempProb))=0;%Set AlphaData to not plot NaNs
% imagesc(-60:60/divider,Range.rkm,tempProb,'AlphaData',imAlpha)
% hold on
% plot(meanTemp,Range.rkm,'b','linewidth',2)
% plot(meanTemp+stdTemp,Range.rkm,'--k')
% plot(meanTemp-stdTemp,Range.rkm,'--k')
% hold off
% colormap(flipud(hot))%colormap hot
% set(gca,'Color','#D3D3D3')
% a = colorbar;
% a.Label.String = 'Occurrences';
% hold on
% xline(0)
% hold off
% ylim([0 4])
% xlim([-30 30])
% set(gca, 'YDir','normal')
% xlabel('\Delta WV % MPD-Sonde (^oC)')
% ylabel('Range (km)')
% title('WV corrected')
% grid on

% %===poisson thin filter in time
% figure(2828)
% title('time')
% subplot(1,6,1)
% semilogx(Counts.Poissonthin.timeWidthon*Time.t_step,Range.rkm,'.')
% xlim([1 10^4])
% grid on
% title('on')
% ylabel('Range')
% subplot(1,6,2)
% semilogx(Counts.Poissonthin.timeWidthoff*Time.t_step,Range.rkm,'.')
% xlim([1 10^4])
% grid on
% title('off')
% subplot(1,6,3)
% semilogx(Counts.Poissonthin.timeWidthon_mol*Time.t_step,Range.rkm,'.')
% xlim([1 10^4])
% grid on
% title('on mol')
% subplot(1,6,4)
% semilogx(Counts.Poissonthin.timeWidthoff_mol*Time.t_step,Range.rkm,'.')
% xlim([1 10^4])
% grid on
% title('of mol')
% subplot(1,6,5)
% semilogx(Counts.Poissonthin.timeWidthwvon*Time.t_step,Range.rkm,'.')
% xlim([1 10^4])
% grid on
% title('on wv')
% subplot(1,6,6)
% semilogx(Counts.Poissonthin.timeWidthwvoff*Time.t_step,Range.rkm,'.')
% xlim([1 10^4])
% grid on
% title('off wv')
% 
% %==poisson thinning filter in range
% figure(2829)
% title('time')
% subplot(6,1,1)
% semilogy(Time.thr,Counts.Poissonthin.rangeWidthon*Range.rangeBin,'.')
% ylim([1 300])
% xlim([Time.thr(1) Time.thr(end)])
% grid on
% title('on')
% subplot(6,1,2)
% semilogy(Time.thr,Counts.Poissonthin.rangeWidthoff*Range.rangeBin,'.')
% ylim([1 300])
% xlim([Time.thr(1) Time.thr(end)])
% grid on
% title('off')
% subplot(6,1,3)
% semilogy(Time.thr,Counts.Poissonthin.rangeWidthon_mol*Range.rangeBin,'.')
% ylim([1 300])
% xlim([Time.thr(1) Time.thr(end)])
% grid on
% title('on_mol')
% subplot(6,1,4)
% semilogy(Time.thr,Counts.Poissonthin.rangeWidthoff_mol*Range.rangeBin,'.')
% ylim([1 300])
% xlim([Time.thr(1) Time.thr(end)])
% grid on
% title('of_mol')
% subplot(6,1,5)
% semilogy(Time.thr,Counts.Poissonthin.rangeWidthwvon*Range.rangeBin,'.')
% ylim([1 300])
% xlim([Time.thr(1) Time.thr(end)])
% grid on
% title('on wv')
% subplot(6,1,6)
% semilogy(Time.thr,Counts.Poissonthin.rangeWidthwvoff*Range.rangeBin,'.')
% ylim([1 300])
% xlim([Time.thr(1) Time.thr(end)])
% grid on
% title('off wv')
% xlabel('Time UTC')

%==sonde and dial difference
figure(8797861)
tempComparison = ones(Range.i_range,size(Sonde.sonde_ind,2));
for jj = 1:size(Sonde.sonde_ind,2)
    tempComparison(:,jj) = diag(Temperature.T_finalm(:,Sonde.sonde_ind(:,jj)))-Sonde.T_sonde(:,jj);
end
if ~isempty(Sonde.sonde_ind)
 meanTemp = mean(tempComparison,2,'omitnan');
else
    meanTemp = nan;
end

plot(meanTemp,Range.rkm,'linewidth',3)
hold on
for iii=1:size(Sonde.sonde_ind,2)
    plot(diag(Temperature.T_finalm(:,Sonde.sonde_ind(:,iii)))-Sonde.T_sonde(:,iii),Range.rkm)
%hold on
end
hold off
grid on
xlim([-10 10])
title('T_{MPD}-T_{sonde}')
xlabel('\Delta T (^oC)')
ylabel('Range (km)')
xline(0)
legend('Mean')


%==sonde and dial difference
figure(8797862)
tempComparison = ones(Range.i_range,size(Sonde.sonde_ind,2));
for jj = 1:size(Sonde.sonde_ind,2)
    tempComparison(:,jj) = (diag(Alpha.alpha_totalm(:,Sonde.sonde_ind(:,jj)))-Sonde.absorption_sonde{jj})./Sonde.absorption_sonde{jj}*100;
end
if ~isempty(Sonde.sonde_ind)
 meanTemp = mean(tempComparison,2,'omitnan');
else
    meanTemp = nan;
end
plot(meanTemp,Range.rkm,'linewidth',3)
hold on
for iii=1:size(Sonde.sonde_ind,2)
    %plot(diag(Temperature.T_finalm(:,Sonde.sonde_ind(:,iii)))-Sonde.T_sonde(:,iii),Range.rkm)
    plot((diag(Alpha.alpha_totalm(:,Sonde.sonde_ind(:,iii)))-Sonde.absorption_sonde{iii})./Sonde.absorption_sonde{iii}*100,Range.rkm)
%hold on
end
hold off
grid on
xlim([-10 10])
title('a_{MPD}-a_{sonde}')
xlabel('\Delta a % (^oC)')
ylabel('Range (km)')
xline(0)
legend('Mean')

% figure(8797862)
% %for iii=1:size(Sonde.sonde_ind,2)
%     for iii=1:5
%     plot(diff(Sonde.T_sonde(:,iii))./Range.rangeBin,Range.rm(1:end-1))
% hold on
%  %%plot(diff(diag(Temperature.T_finalm(:,Sonde.sonde_ind(:,iii))))./Range.rangeBin,Range.rm(1:end-1),'--')
% end
% xline(-.0065)
% hold off
% grid on
% xlim([-.05 .05])


figure
errorbar(diag(Temperature.T_finalm(:,p_point))-273,Range.rkm,diag(Temperature.TempStds(:,p_point)))
errorbar(Range.rkm,diag(Temperature.T_finalm(:,p_point))-273,diag(Temperature.TempStds(:,p_point)))
hold on
plot(Range.rkm,Sonde.T_sonde(:,sonde_index)-273)
xlabel('Range (km)')
ylabel('Temperature (^oC)')
legend('T_{MPD}','T_{radiosonde}')
title(sprintf(['Temperature\n' datestr(Time.date_ts(p_point(1)))]))
xlim([0 4])
view([90 -90])
grid on

% for ii=1:size(Sonde.T_sonde,2)
%     p_point2 = Sonde.sonde_ind(:,ii);
%     figure
%     errorbar(Range.rkm,diag(Temperature.T_finalm(:,p_point2))-273,diag(Temperature.TempStds(:,p_point2)))
%     hold on
%     plot(Range.rkm,Sonde.T_sonde(:,ii)-273)
%     xlabel('Range (km)')
%     ylabel('Temperature (^oC)')
%     legend('T_{MPD}','T_{radiosonde}')
%     title(sprintf(['Temperature\n' datestr(Time.date_ts(p_point2(1)))]))
%     xlim([0 4])
%     view([90 -90])
%     grid on
%     saveas(gcf,['C:\Users\Owen\OneDrive - Montana State University\Research\Reports\11_22_22\temp all sondes\' sprintf('temp%d.png',ii)]);
% end

figure(1114)
subplot(1,2,1)
T_pot = Temperature.T_finalm.*(0.986923./Model.P).^0.286;
T_pot_sonde = Sonde.T_sonde.*(0.986923./Sonde.P_sonde).^0.286;
plot(T_pot(:,p_point(1)),Range.rkm)
hold on
plot(T_pot_sonde(:,sonde_index),Range.rkm)
hold off
xlabel('Potential temperature')
ylabel('Range')
grid on
subplot(1,2,2)
T_pot_diff = diff(T_pot,1);
T_pot_sdiff = diff(T_pot_sonde,1);
plot(T_pot_diff(:,p_point(1))./Range.rangeBin,Range.rkm(1:end-1)+Range.rangeBin/2/1000)
hold on
plot(T_pot_sdiff(:,sonde_index)./Range.rangeBin,Range.rkm(1:end-1)+Range.rangeBin/2/1000)
hold off
xlabel('diff Potential temperature')
grid on

figure(11115)
subplot(2,1,1)
h = pcolor(Time.date_ts,Range.rkm,T_pot);
set(h, 'EdgeColor', 'none');
subplot(2,1,2)
h = pcolor(Time.date_ts,Range.rkm(1:end-1),T_pot_diff);
set(h, 'EdgeColor', 'none');



% figure(1234)
% plot(Time.date_ts(1:12082),Model.Ts(1:12082)-273,'k','linewidth',2)
% hold on
% plot(Time.date_ts(12083:end),Model.Ts(12083:end)-273,'k','linewidth',2)
% grid on
% ylabel('Suface Temperature (^oC)')




end

function timeAxis(Format,Time)
%xlabel('Time UTC')
xticks(datenum(Time.date_ts(1)): Format.tickSpacing :datenum(Time.date_ts(end)))
%xticks(datenum(datetime(2022,04,30,12,0,0,'TimeZone','UTC')): Format.tickSpacing :datenum(Time.date_ts(end)))
%datenum(datetime(2022,04,30,14,0,0,'TimeZone','UTC'))
xtickangle(Format.tickAngle)
datetick('x',Format.dateTickFormat,'keeplimits','keepticks')


%xtickformat('MM-dd')
end

