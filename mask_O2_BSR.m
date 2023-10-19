function [SNRm , cloud_SDm_above, cloud_SDm,o2on_SNR] = mask_O2_BSR(cloud_p_point,SNR_threshold,SD_threshold,~,~,Counts,BGmult,Time,Range,BSR,lowAlt)
%File: mask_O2.m
%Date: 03/16/2020
%Author: Owen Cruikshank
%Inputs:
%   -o2on:[photons] (range x time) vector of online return photon counts
%   -o2off:[photons] (range x time) vector of offline return photon counts
%   -rm:[m] (range x 1) vector of range dimension
%   -ts:[s] (1 x time) vector of time dimension
%   -cloud_p_point:[hr] point to plot cloud plots. If 0 will not plot.
%   -SNR_threshold: Signal to noise rato threshold for SNR mask
%   -SD_threshold: Standard deviation threshold for cloud mask
%   -oversample: Number of range bin averaging
%   -t_avg: Number of time bin averaging
%
%Outputs:
%   -SNRm:[none] (range x time) calculated SNR mask
%   -cloud_SDm_above:[none] (range x time) calculated cloud mask
%   -cloud_SDm:[none] (range x time) calculated cloud mask with only clouds
%   -o2on_SNR:[none] (range x time) signal to noise ratio for online counts


% ===================
% === SNR O2 mask ===
% ===================

o2on_SNR = sqrt(Counts.NBins).*Counts.o2on./(sqrt(Counts.o2on+Counts.bg_o2on*BGmult));   % Calculate shot noise SNR with background and with averages
o2on_SNRm = ones(size(o2on_SNR));       % Initialize matrix
o2on_SNRm(o2on_SNR < SNR_threshold & Range.rm > 1000 | isnan(Counts.o2off) | (Range.rm.*ones(size(o2on_SNR)))<lowAlt) = 0;% Replace low SNR points with zeros

o2off_SNR = sqrt(Counts.NBins).*Counts.o2off./(sqrt(Counts.o2off+Counts.bg_o2off*BGmult));   % Calculate shot noise SNR with background and with averages
o2off_SNRm = ones(size(o2off_SNR));
o2off_SNRm((o2off_SNR < SNR_threshold & Range.rm > 1000) | isnan(Counts.o2off) | (Range.rm.*ones(size(o2off_SNR)))<lowAlt) = 0;

SNRm = o2off_SNRm .* o2on_SNRm;         %Final combination SNR matrix mask with 0 in places of low snr


% ===================
% === Cloud O2 mask ===
% ===================

%Range corrected returns
% o2on_range_corrected = Counts.o2on.*Range.rm.^2;
% o2off_range_corrected = Counts.o2off.*Range.rm.^2;


%standard deviation
%stdNeiborhood = true(7,5);   
stdNeiborhood = true(3,1); % NxN window arround each point

% o2on_SD = stdfilt(o2on_range_corrected,stdNeiborhood);  % Use a square standard deviation filter on the range corrected return to find a standard deviation
% cloud_SDm_on = ones(size(Counts.o2on));                        % Initialize mask matrix
% cloud_SDm_on(o2on_SD > SD_threshold) = 0;                       % Set points with a standard deviation greater than the threshold to 0 inicating a cloudy area
% 
% o2off_SD = stdfilt(o2off_range_corrected,stdNeiborhood);
% cloud_SDm_off = ones(size(Counts.o2on));
% cloud_SDm_off(o2off_SD > SD_threshold) = 0;
% 
% cloud_SDm = cloud_SDm_off .* cloud_SDm_on;              % Combine on and offline cloud masks


%USE BSR as cloud mask
BSR_SD = stdfilt(BSR,stdNeiborhood);  % Use a square standard deviation filter on the range corrected return to find a standard deviation
cloud_SDm_BSR = ones(size(BSR));                        % Initialize mask matrix
%cloud_SDm_BSR(BSR_SD > SD_threshold & Range.rm > 500 | isnan(Counts.o2off)) = 0;                       % Set points with a standard deviation greater than the threshold to 0 inicating a cloudy area

cloud_SDm_BSR((BSR_SD > SD_threshold) & (Range.rm > lowAlt) | isnan(Counts.o2off)) = 0;  

cloud_SDm = cloud_SDm_BSR;


%Set values above the cloud to -1
cloud_SDm_above = cloud_SDm;
diff_cloud_SDm = diff(cloud_SDm);                           % Calculate difference along range axis
diff_SNRm = diff(SNRm);

[cloud_index_row,cloud_index_col] = find(diff_cloud_SDm);   % Find non-zero values in difference
if ~isempty(cloud_index_row)                                 % Check if any clouds exist
    cloud_index = zeros(1,length(Time.ts));                          % Initalize cloud range index matrix

    index_matrix = ones(length(Range.rm),length(Time.ts)) .* (1:length(Range.rm))'; % Create matrix to match range indecies

    % Set range index value to index of cloud
    cloud_index(cloud_index_col(1))=cloud_index_row(1);         
    for i = 2:length(cloud_index_col)
        if cloud_index_col(i)~=cloud_index_col(i-1)
            cloud_index(cloud_index_col(i))=cloud_index_row(i);
        end
    end
    
    cloud_SDm_above(index_matrix > cloud_index & cloud_index ~= 0)=-1; % Set all values above cloud index to -1
    cloud_SDm_above = cloud_SDm_above .* cloud_SDm;                    % Create final matrix where clouds are represended by zeros and all data above represented by -1   
end

[SNR_index_row,SNR_index_col] = find(diff_SNRm);   % Find non-zero values in difference
if ~isempty(SNR_index_row)                                 % Check if any clouds exist
    % Initalize cloud range index matrix
    SNR_index = ones(1,length(Time.ts))*length(Range.rm);
    
    %index_matrix = ones(length(Range.rm),length(Time.ts)) .* (1:length(Range.rm))'; % Create matrix to match range indecies

    % Set range index value to index of cloud
    SNR_index(SNR_index_col(1))=SNR_index_row(1);         
    for i = 2:length(SNR_index_col)
        if SNR_index_col(i)~=SNR_index_col(i-1)
            SNR_index(SNR_index_col(i))=SNR_index_row(i);
        end
    end
end

%%

%---Plotting
if cloud_p_point ~= 0
[~,p_point] = min(abs(cloud_p_point-Time.ts/60/60));%find closest value to 338min for comparison to other program
thr = Time.ts/60/60;
 
figure(565)
subplot(2,1,1)
%plot(Range.rm,o2off_SD(:,p_point))
%hold on
%plot(Range.rm,o2on_SD(:,p_point))
plot(Range.rm,BSR_SD(:,p_point),'--')
yline(SD_threshold)
%hold off
legend('on','off')
title('Standard Deviation filter')
xlabel('Range(m)')

subplot(2,1,2)
plot(Range.rm,o2on_SNR(:,p_point))
hold on
plot(Range.rm,o2off_SNR(:,p_point))
hold off
yline(SNR_threshold)
legend('on','off')
title('Signal to noise ratio')
xlabel('Range(m)')


figure(60)
subplot(3,1,1)
imagesc(thr,Range.rm,o2on_SNRm)
hold on
xline(thr(p_point),'r');
hold off
set(gca, 'YDir','normal')
title('o2on_SNRm')
subplot(3,1,2)
imagesc(thr,Range.rm,o2off_SNRm)
hold on
xline(thr(p_point),'r');
hold off
set(gca, 'YDir','normal')
title('oo2off_SNRm')

subplot(3,1,3)
imagesc(thr,Range.rm,SNRm)
hold on
xline(thr(p_point),'r');
hold off
set(gca, 'YDir','normal')
title('SNRm')

figure(61)
imagesc(thr,Range.rm,cloud_SDm)
hold on
xline(thr(p_point),'r');
hold off
set(gca, 'YDir','normal')
title('cloudSDM')

end
end