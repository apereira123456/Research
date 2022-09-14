%% Clear
close all; clear all; clc;

%% Load Library
load('phaseLibrary.mat')

%% Setup Parameters
% Specify excitation laser wavelength in nm
lambda_excitation = 532;

first_guess = struct('b4c1', [539.63 1.55 541.33 1.28], 'b4c2', [546.00 0.32 547.52 0.37], ...
                     'b4c3', [553.33 1.55 555.58 1.55 556.68 0.96 557.50 0.64], ...
                     'b4c4', [560.79 0.80 561.98 1.28 564.66 2.51], 'si', [547.4 1], ...
                     'sic', [560.76 0.283]);

% Specify spectral bands to be analyzed
wave_start = 535.74;
wave_end = 595;

%% Data Selection Prompt
[file, path] = uigetfile({'*.txt'});

data_name = fullfile(path, file);
data = table2array(readtable(data_name));

%% Determine data array size
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

%% Phase ID
% Score Map
scoreMap(:,:,1) = sidsam(dcube, phaseLibrary(1).Reflectance);
scoreMap(:,:,2) = sidsam(dcube, phaseLibrary(2).Reflectance);
scoreMap(:,:,3) = sidsam(dcube, phaseLibrary(3).Reflectance);
scoreMap(:,:,4) = sidsam(dcube, phaseLibrary(4).Reflectance);

% Class Map
[~,classMap] = min(scoreMap,[],3);

figure
imagesc(classMap)

ticks = linspace(1,4,4);
colormap(jet(4))
colorbar('Ticks',ticks,'TickLabels',[phaseLibrary.Name])
clim([0.5 4.5])

%% Peak Fit
fit_data = cell(y_size,x_size);

for j = 1:y_size
    for i = 1:x_size
        if classMap(j,i) == 1
            [~,b4c1_start] = min(abs(wavelength - 536.3));
            [~,b4c1_end] = min(abs(wavelength - 544.0));
            
            [~,b4c2_start] = min(abs(wavelength - 545.0));
            [~,b4c2_end] = min(abs(wavelength - 548.8));
            
            [~,b4c4_start] = min(abs(wavelength - 558.9));
            [~,b4c4_end] = min(abs(wavelength - 571.53));
            
            signal1 = [wavelength(b4c1_start:b4c1_end) squeeze(dcube(j,i,b4c1_start:b4c1_end))];
            signal2 = [wavelength(b4c2_start:b4c2_end) squeeze(dcube(j,i,b4c2_start:b4c2_end))];
            signal4 = [wavelength(b4c4_start:b4c4_end) squeeze(dcube(j,i,b4c4_start:b4c4_end))];
            
            [FitResults1, ~] = peakfit(signal1,0,0,2,1,0,3,first_guess.b4c1,1,0,0);         
            [FitResults2, ~] = peakfit(signal2,0,0,2,1,0,3,first_guess.b4c2,1,0,0);         
            [FitResults4, ~] = peakfit(signal4,0,0,3,1,0,3,first_guess.b4c4,1,0,0);

            fit_data(j,i) = mat2cell([FitResults1; FitResults2; FitResults4], 7, 5);
        
        elseif classMap(j,i) == 4
            [~,sic_start] = min(abs(wavelength - 557.97));
            [~,sic_end] = min(abs(wavelength - 563.57));
            sic_signal = [wavelength(sic_start:sic_end) squeeze(dcube(j,i,sic_start:sic_end))];

            [FitResults, ~] = peakfit(sic_signal,0,0,1,1,0,3,first_guess.sic,1,0,0);            
            fit_data(j,i) = mat2cell(FitResults, 1, 5);
        end
    end
end

%%
save('fit_data')

%%
load('fit_data.mat')

%% Stress Map
peakMap = zeros(y_size,x_size);
for j = 1:y_size
    for i = 1:x_size
        if classMap(j,i) == 1
            [row, ~] = find(fit_data{j,i}(:,1) == 3);
            peakMap(j,i) = fit_data{j,i}(row,2);

%             peakMap(j,i) = fit_data{j,i}(3,2);
        else
            peakMap(j,i) = NaN;
        end
    end
end

% for j = 1:y_size
%     for i = 1:x_size
%         if classMap(j,i) == 4
%             peakMap(j,i) = fit_data{j,i}(1,2);
%         else
%             peakMap(j,i) = NaN;
%         end
%     end
% end

figure
heatmap = imagesc(peakMap);
set(heatmap,'AlphaData',~isnan(peakMap))
colormap("jet")
colorbar
clim([564.7 564.8])
% clim([545.945 545.985])
% clim([560.87 560.92])

%%
% Generate the new figure name
figure_name = append(path, 'Raman Stress Map.png');

% Save the figure
export_fig(figure_name,'-transparent','-m5');