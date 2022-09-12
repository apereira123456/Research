%% Clear
close all; clear all; clc;

%% Load Library
load('phaseLibrary.mat')

%% Setup Parameters
% Specify excitation laser wavelength in nm
lambda_excitation = 532;

% first_guess = struct('b4c', [270 60 320 60 480 20 530 15 724 50 825 50 1000 40 1089 90], ...
%                      'sic1', [554.60 0.170 554.88 0.227 555.31 0.227 555.52 0.227], );

first_guess = struct('sic', [560.76 0.283]);

% Specify spectral bands to be analyzed
wave_start = 535.74;
wave_end = 595;

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
% Store Raman shift
wavelength = data(s_size:-1:1,3);

% Convert Raman shift to Raman wavelength
wavelength = ( (lambda_excitation)^(-1) - (wavelength ./ 10^7) ).^(-1);

% Initialize the cell that stores the Raman data
dcube = zeros(y_size,x_size,s_size);

for j = 1:y_size
    for i = 1:x_size
        a = (i-1)*s_size + (j-1)*x_size*s_size + 1;
        b = i*s_size + (j-1)*x_size*s_size;

        dcube(j,i,:) = data(b:-1:a,4);
    end
end

[~, start_index] = min(abs(wavelength - wave_start));
[~, end_index] = min(abs(wavelength - wave_end));

wavelength = wavelength(start_index:end_index);
dcube = normalize(dcube(:,:,start_index:end_index),3,'range');

%% Viewer
% hcube = hypercube(dcube, wavelength);
% hyperspectralViewer(hcube);

%% Peak Fit
% SiC
[~,peak_start] = min(abs(wavelength - 557.97));
[~,peak_end] = min(abs(wavelength - 563.57));

fit_data = cell(y_size,x_size);

for j = 1:y_size
    for i = 1:x_size
        signal = [wavelength(peak_start:peak_end) squeeze(dcube(j,i,peak_start:peak_end))];
        [FitResults,FitError] = peakfit(signal,0,0,1,1,0,3,first_guess.sic,1,0,0);
        
        fit_data(j,i) = mat2cell(FitResults, 1, 5);
    end
end

%%
save('fit_data')

%%
% load('fit_data.mat')

%% Stress Map
peakMap = zeros(y_size,x_size);
for j = 1:y_size
    for i = 1:x_size
        peakMap(j,i) = fit_data{j,i}(1,2);
    end
end

figure
heatmap = imagesc(peakMap);
set(heatmap,'AlphaData',~isnan(peakMap))
colormap("jet")
colorbar
% clim([0 0.11])

%%
% Generate the new figure name
figure_name = append(path, 'Raman Stress Map.png');

% Save the figure
export_fig(figure_name,'-transparent','-m5');