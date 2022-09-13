%% Clear
close all; clear all; clc;

%% Load Library
load('phaseLibrary.mat')

%% Setup Parameters
% Specify excitation laser wavelength in nm
lambda_excitation = 532;

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

hcube = hypercube(dcube, wavelength);

%% Viewer
% hyperspectralViewer(hcube);

%% Match Spectra
%% SIDSAM
scoreMap(:,:,1) = sidsam(dcube, phaseLibrary(1).Reflectance);
scoreMap(:,:,2) = sidsam(dcube, phaseLibrary(2).Reflectance);
scoreMap(:,:,3) = sidsam(dcube, phaseLibrary(3).Reflectance);
scoreMap(:,:,4) = sidsam(dcube, phaseLibrary(4).Reflectance);

% %% B4C
% figure
% title('Boron Carbide')
% imagesc(scoreMap(:,:,1))
% colormap(jet);
% colorbar
% clim([0 150])
% 
% %% Si
% figure
% title('Silicon')
% imagesc(scoreMap(:,:,2))
% colormap(jet);
% colorbar
% clim([0 3500])
% 
% %% SiB6
% figure
% title('Silicon Boride')
% imagesc(scoreMap(:,:,3))
% colormap(jet);
% colorbar
% clim([0 200])
% 
%% SiC
figure
title('SiC')
imagesc(scoreMap(:,:,4))
colormap(jet);
colorbar
clim([0 0.7])

% scoreMap(:,:,1) = jmsam(dcube, phaseLibrary(1).Reflectance);
% scoreMap(:,:,2) = jmsam(dcube, phaseLibrary(2).Reflectance);
% scoreMap(:,:,3) = jmsam(dcube, phaseLibrary(3).Reflectance);
% scoreMap(:,:,4) = jmsam(dcube, phaseLibrary(4).Reflectance);

%% Class Map
[~,classMap] = min(scoreMap,[],3);

figure
imagesc(classMap)

ticks = linspace(1,4,4);
colormap(jet(4))
colorbar('Ticks',ticks,'TickLabels',[phaseLibrary.Name])
clim([0.5 4.5])