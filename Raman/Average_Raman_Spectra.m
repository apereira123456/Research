%% Clear
close all; clear all; clc;

%% Setup Parameters
% Specify excitation laser wavelength in nm
lambda_excitation = 532;

% Specify spectral bands to be analyzed
wave_start = 130;                   % spectral bands to be analyzed
wave_end = 1500;                    % spectral bands to be analyzed

%% Data Selection Prompt
[file, path] = uigetfile({'*.txt'});

UIFigure.Visible = 'on';
DataEditField.Value = file;

data_name = fullfile(path, file);
data = readmatrix(data_name);


%% Determine data array size
s_size = find(data(2:end,3) == data(1,3), 1, 'first');

x_size = round((data(end,1) - data(1,1)) / (data(s_size + 1,1) - data(1,1)) + 1);
y_size = round((data(end,2) - data(1,2)) / (data(x_size*s_size + 1,2) - data(1,2)) + 1);

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

%%
avg = flipud(squeeze(mean(dcube, [1 2])));
SiBC = normalize(avg, 'range');

%%
Wavenumber = flipud(wavenumber);
T = table(Wavenumber, SiBC);
writetable(T, 'Standard Raman Spectra 1800.txt', 'Delimiter', 'tab')
