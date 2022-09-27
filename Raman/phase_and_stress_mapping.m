%% Clear
close all; clear all; clc

lambda_excitation = 532;            % excitation laser wavelength in nm

wave_start = 130;                   % spectral bands to be analyzed
wave_end = 1500;                    % spectral bands to be analyzed

first_guess = struct('b4c1', [265.8 54.5 324.1 45.1], 'b4c2', [482.0 11.3 532.7 13.2], ...
                     'b4c3', [724.4 54.5 797.7 54.5 833.5 33.8 859.8 22.6], ...
                     'b4c4', [965.0 28.2 1002.6 45.1 1087.2 88.3], 'si', [528.8 35.3], ...
                     'sic', [964.1 10.0]);

peak_window = struct('b4c1', [150.7 414.6], 'b4c2', [448.4 575.4], ...
                     'b4c4', [904.7 1300.1], 'si', [448.4 615.2], ...
                     'sic', [874.9 1053.0]);

%% Standard
standard = readtable('Standard Raman Spectra.txt');

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

for j = 1:y_size
    for i = 1:x_size
        a = (i-1)*s_size + (j-1)*x_size*s_size + 1;
        b = i*s_size + (j-1)*x_size*s_size;

        dcube(j,i,:) = data(b:-1:a,4);
    end
end

[~, start_index] = min(abs(wavenumber - wave_start));
[~, end_index] = min(abs(wavenumber - wave_end));

num_bands = end_index - start_index + 1;

wavenumber = wavenumber(start_index:end_index);
dcube = normalize(dcube(:,:,start_index:end_index),3,'range');


%% Match Spectra
score_map(:,:,1) = sidsam(dcube, standard.B4C);
score_map(:,:,2) = sidsam(dcube, standard.Si);
score_map(:,:,3) = sidsam(dcube, standard.SiB6);
score_map(:,:,4) = sidsam(dcube, standard.SiC);

%% Phase Map
[~,phase_map] = min(score_map,[],3);

figure
imagesc(phase_map)

ticks = linspace(1,4,4);
labels = standard.Properties.VariableNames;

colormap jet(4);
colorbar('Ticks',ticks, 'TickLabels',labels(2:end));
clim([0.5 4.5])

pause(0.5)


%% Parallel Processing
if isempty(gcp('nocreate'))
    parpool('local')
end


%% Fit Peaks
peak_fit = cell(y_size,x_size);

% Add code to remove peak fits with poor goodness of fit

tic
for j = 1:y_size
    fprintf('%d\n',j);
    parfor i = 1:x_size
        if phase_map(j,i) == 1          
            [~,b4c1_start] = min(abs(wavenumber - peak_window.b4c1(1)));
            [~,b4c1_end] = min(abs(wavenumber - peak_window.b4c1(2)));
            
            [~,b4c2_start] = min(abs(wavenumber - peak_window.b4c2(1)));
            [~,b4c2_end] = min(abs(wavenumber - peak_window.b4c2(2)));
            
            [~,b4c4_start] = min(abs(wavenumber - peak_window.b4c4(1)));
            [~,b4c4_end] = min(abs(wavenumber - peak_window.b4c4(2)));
            
            signal1 = [wavenumber(b4c1_start:b4c1_end) squeeze(dcube(j,i,b4c1_start:b4c1_end))];
            signal2 = [wavenumber(b4c2_start:b4c2_end) squeeze(dcube(j,i,b4c2_start:b4c2_end))];
            signal4 = [wavenumber(b4c4_start:b4c4_end) squeeze(dcube(j,i,b4c4_start:b4c4_end))];
            
            [FitResults1, ~] = peakfit(signal1,0,0,2,1,0,3,first_guess.b4c1,1,0,0);         
            [FitResults2, ~] = peakfit(signal2,0,0,2,1,0,3,first_guess.b4c2,1,0,0);         
            [FitResults4, ~] = peakfit(signal4,0,0,3,1,0,3,first_guess.b4c4,1,0,0);

            fit_data(j,i) = mat2cell([FitResults1(:,2:5); FitResults2(:,2:5); FitResults4(:,2:5)], 7, 4);         
        
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

peak_fit = fit_data;


%% Preallocate Vectors
total(:,:,1) = sum(8 * (phase_map == 1),'all');
total(:,:,2) = sum(2 * (phase_map == 2),'all');
total(:,:,3) = sum(2 * (phase_map == 3),'all');
total(:,:,4) = sum(2 * (phase_map == 4),'all');
total(:,:,5) = sum(2 * (phase_map == 0),'all');
num_rows = sum(total,'all');

[M, X, Y, H, Position, Height, FWHM, Area] = deal(NaN(num_rows,1));

table_row = 1;

for j = 1:y_size
    for i = 1:x_size
        if phase_map(j,i) == 1         
            M(table_row) = 1;
            X(table_row) = i;
            Y(table_row) = j;
            H(table_row) = 7;
            Position(table_row:table_row + 6) = peak_fit{j,i}(:,1);
            Height(table_row:table_row + 6) = peak_fit{j,i}(:,2);
            FWHM(table_row:table_row + 6) = peak_fit{j,i}(:,3);
            Area(table_row:table_row + 6) = peak_fit{j,i}(:,4);

            table_row = table_row + 8;           
        
        elseif phase_map(j,i) == 4
            M(table_row) = 4;
            X(table_row) = i;
            Y(table_row) = j;
            H(table_row) = 1;
            Position(table_row) = peak_fit{j,i}(:,1);
            Height(table_row) = peak_fit{j,i}(:,2);
            FWHM(table_row) = peak_fit{j,i}(:,3);
            Area(table_row) = peak_fit{j,i}(:,4);

            table_row = table_row + 2;
            
        else
            M(table_row) = 0;
            X(table_row) = i;
            Y(table_row) = j;
            H(table_row) = 0;

            table_row = table_row + 2;
        end
    end
end


%% Write Table
peak_export = table(M, X, Y, H, Position, Height, FWHM, Area);
writetable(peak_export, 'peak_export.txt', 'Delimiter', 'tab')

%% Save
save('fit_data')

%% Load
load('fit_data.mat')