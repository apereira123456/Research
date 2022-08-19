%% Clear
close all; clear all; clc;

b4c_standard = table2array(readtable('Standard Spectra/B4C.txt'));
si_standard = table2array(readtable('Standard Spectra/Si.txt'));
sib6_standard = table2array(readtable('Standard Spectra/SiB6.txt'));

standard(:,1) = b4c_standard(:,3);
standard(:,2) = b4c_standard(:,4);
standard(:,3) = si_standard(:,4);
standard(:,4) = sib6_standard(:,4);

% Determine the indices associated with the truncation start and stop
trunc_range = [150, 1500];
[~,trunc_start] = min(abs(standard(:,1) - trunc_range(1)));
[~,trunc_end] = min(abs(standard(:,1) - trunc_range(2)));

standard = standard(trunc_start: -1: trunc_end,:);

%% Data Selection Prompt
[file, path] = uigetfile({'*.txt'});

data_name = fullfile(path, file);
data = table2array(readtable(data_name));

%% Logic for Scans of Different Dimensions
% Determine data array size
[~,x,~] = unique(data(:,1),'rows');
[~,y,~] = unique(data(:,2),'rows');
[~,s,~] = unique(data(:,3),'rows');

[x_size, ~] = size(x);
[y_size, ~] = size(y);
[s_size, ~] = size(s);

%% Convert raw data to usable format
% Initialize the cell that stores the Raman data
map_data = cell(x_size,y_size,2);

% Determine the indices associated with the truncation start and stop
[~,trunc_start] = min(abs(data(1:s_size,3) - trunc_range(2)));
[~,trunc_end] = min(abs(data(1:s_size,3) - trunc_range(1)));

trunc_size = trunc_end - trunc_start + 1;

% Adjust the truncation start and stop so the values can be used with indexing variables a and b
start_offset = trunc_start - 1;
end_offset = s_size - trunc_end;

for j = 1:y_size
    for i = 1:x_size
        a = (i-1)*s_size + (j-1)*x_size*s_size + 1;
        a = a + start_offset;

        b = i*s_size + (j-1)*x_size*s_size;
        b = b - end_offset;

        % Store the truncated data and normalize the intensities
        map_data(i,j,1) = mat2cell(data(b:-1:a,3),trunc_size,1);
        map_data(i,j,2) = mat2cell(normalize(data(b:-1:a,4),1,'range'),trunc_size,1);
    end
end

% for j = 1:y_size
%     for i = 1:x_size
%         
%     end
% end

% [b4c_corr,b4c_lag] = xcorr(map_data{1,1,2}(1:end-1,1),standard(:,2));
% [si_corr,si_lag] = xcorr(map_data{1,1,2}(1:end-1,1),standard(:,3));
% [sib6_corr,sib6_lag] = xcorr(map_data{1,1,2}(1:end-1,1),standard(:,4));


wav = getMatchedFilter(waveform);
filter = phased.MatchedFilter('Coefficients',wav);
taylorfilter = phased.MatchedFilter('Coefficients',wav,'SpectrumWindow','Taylor');

[P1,Q1] = rat(map_data{55,24,2}/standard(:,1)Fs/Fs1);          % Rational fraction approximation
T1 = resample(T1,P1,Q1);        % Change sampling rate by rational factor
T2 = resample(T2,P2,Q2);        % Change sampling rate by rational factor

[b4c_corr,b4c_lag] = xcorr(map_data{55,24,2}(1:end-1,1),standard(:,2));
[si_corr,si_lag] = xcorr(map_data{55,24,2}(1:end-1,1),standard(:,3));
[sib6_corr,sib6_lag] = xcorr(map_data{55,24,2}(1:end-1,1),standard(:,4));

figure
ax(1) = subplot(3,1,1); 
plot(b4c_lag,b4c_corr)
ylabel('Intensity')
grid on
title('Cross-correlation between Map and B4C Standard')

ax(2) = subplot(3,1,2); 
plot(si_lag,si_corr)
ylabel('Intensity') 
grid on
title('Cross-correlation between Map and Si Standard')

ax(3) = subplot(3,1,3); 
plot(sib6_lag,sib6_corr)
ylabel('Intensity')
grid on
title('Cross-correlation between Map and SiB6 Standard')

xlabel('Wavenumber (cm^-1)') 
% axis(ax(1:3),[150 1500 -1 1 ])