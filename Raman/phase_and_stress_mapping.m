%% Clear
close all; clear all; clc;

delete(gcp('nocreate'));
pool = parpool('local');

%% Setup Parameters
% Specify excitation laser wavelength in nm
lambda_excitation = 532;

first_guess = struct('b4c1', [482.0 11.3 532.7 13.2], 'b4c2', [1087.2 88.3], ...
                     'si', [528.8 35.3], 'sic', [964.1 10.0], ...
                     'sibc1', [482.0 11.3 532.7 13.2], 'sibc2', [1065 88.3]);

peak_window = struct('b4c1', [448.4 575.4], 'b4c2', [1060 1120], ...
                     'si', [448.4 615.2], 'sic', [874.9 1053.0], ...
                     'sibc1', [448.4 575.4], 'sibc2', [1020 1100]);

% Specify spectral bands to be analyzed
wave_start = 130;                   % spectral bands to be analyzed
wave_end = 1500;                    % spectral bands to be analyzed

%% Standard Selection Prompt
UIFigure.Visible = 'off';
[file, temp_path] = uigetfile({'*.txt'});

UIFigure.Visible = 'on';
StandardEditField.Value = file;

standard_name = fullfile(temp_path, file);
standard = readtable(standard_name);

standard = flipud(standard);

[~, start_index] = min(abs(standard.Wavenumber - wave_start));
[~, end_index] = min(abs(standard.Wavenumber - wave_end));

standard = standard(start_index:end_index, :);     

%% Data Selection Prompt
[file, path] = uigetfile({'*.txt'});

UIFigure.Visible = 'on';
DataEditField.Value = file;

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
wavenumber = data(s_size:-1:1,3);    

% Initialize the cell that stores the Raman data
dcube = zeros(y_size,x_size,s_size);

tic
for j = 1:y_size
    for i = 1:x_size
        a = (i-1)*s_size + (j-1)*x_size*s_size + 1;
        b = i*s_size + (j-1)*x_size*s_size;

        dcube(j,i,:) = data(b:-1:a,4);
    end
end
toc

[~, start_index] = min(abs(wavenumber - wave_start));
[~, end_index] = min(abs(wavenumber - wave_end));

wavenumber = wavenumber(start_index:end_index);
dcube = normalize(dcube(:,:,start_index:end_index),3,'range');

%% Viewer
% hcube = hypercube(dcube, wavelength);
% hyperspectralViewer(hcube);

%% Match Spectra
score_map(:,:,1) = sidsam(dcube, standard.B4C);
score_map(:,:,2) = sidsam(dcube, standard.Si);
score_map(:,:,3) = sidsam(dcube, standard.SiB6);
score_map(:,:,4) = sidsam(dcube, standard.SiC);

%% Phase Map
[~,phase_map] = min(score_map,[],3);

%% Peak Fit
fit_data = cell(y_size,x_size);

tic
for j = 1:y_size
    sprintf('%d', j)
    parfor i = 1:x_size
        if phase_map(j,i) == 1                      
            [~,b4c1_start] = min(abs(wavenumber - peak_window.b4c1(1)));
            [~,b4c1_end] = min(abs(wavenumber - peak_window.b4c1(2)));
            
            [~,b4c4_start] = min(abs(wavenumber - peak_window.b4c4(1)));
            [~,b4c4_end] = min(abs(wavenumber - peak_window.b4c4(2)));
            
            signal1 = [wavenumber(b4c1_start:b4c1_end) squeeze(dcube(j,i,b4c1_start:b4c1_end))];
            signal4 = [wavenumber(b4c4_start:b4c4_end) squeeze(dcube(j,i,b4c4_start:b4c4_end))];
            
            [FitResults1, ~] = peakfit(signal1,0,0,2,1,0,3,first_guess.b4c1,1,0,0);         
            [FitResults4, ~] = peakfit(signal4,0,0,3,1,0,3,first_guess.b4c4,1,0,0);

            fit_data(j,i) = mat2cell([FitResults1(:,2:5); FitResults2(:,2:5)], 3, 4);          
        
        elseif phase_map(j,i) == 4
            [~,sic_start] = min(abs(wavenumber - peak_window.sic(1)));
            [~,sic_end] = min(abs(wavenumber - peak_window.sic(2)));

            sic_signal = [wavenumber(sic_start:sic_end) squeeze(dcube(j,i,sic_start:sic_end))];

            [FitResults, ~] = peakfit(sic_signal,0,0,1,1,0,3,first_guess.sic,1,0,0);

            fit_data(j,i) = mat2cell(FitResults(:,2:5), 1, 4);
        end
    end
end
toc

save('fit_data')

% %%
% load('fit_data.mat')
% 
% %% Stress Map
% peakMap = zeros(y_size,x_size);
% for j = 1:y_size
%     for i = 1:x_size
%         if classMap(j,i) == 1
%             [row, ~] = find(fit_data{j,i}(:,1) == 3);
%             peakMap(j,i) = fit_data{j,i}(row,2);
% 
% %             peakMap(j,i) = fit_data{j,i}(3,2);
%         else
%             peakMap(j,i) = NaN;
%         end
%     end
% end
% 
% % for j = 1:y_size
% %     for i = 1:x_size
% %         if classMap(j,i) == 4
% %             peakMap(j,i) = fit_data{j,i}(1,2);
% %         else
% %             peakMap(j,i) = NaN;
% %         end
% %     end
% % end
% 
% figure
% heatmap = imagesc(peakMap);
% set(heatmap,'AlphaData',~isnan(peakMap))
% colormap("jet")
% colorbar
% clim([564.7 564.8])
% % clim([545.945 545.985])
% % clim([560.87 560.92])
% 
% %%
% % Generate the new figure name
% figure_name = append(path, 'Raman Stress Map.png');
% 
% % Save the figure
% export_fig(figure_name,'-transparent','-m5');