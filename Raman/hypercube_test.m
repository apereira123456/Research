%% Clear
close all; clear all; clc;

%% Load Library
load('phaseLibrary.mat')

%% Setup Parameters
% Specify excitation laser wavelength in nm
lambda_excitation = 532;

% Specify spectral bands to remove
wlrange = [500 535.6; 595 633];

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

        dcube(j,i,:) = normalize(data(b:-1:a,4),1,'range');
    end
end

hcube = hypercube(dcube, wavelength);
hcube = removeBands(hcube,'Wavelength', wlrange);

% hyperspectralViewer(newhcube);

%% Match Spectra
% scoreMap = spectralMatch(phaseLibrary, hcube);
% 
% figure
% montage(scoreMap,'Size',[2 3],'BorderSize',10)
% title('Score Map Obtained for Each Pure Spectrum','FontSize',14)
% colormap(jet);
% colorbar
% clim([0 5])

scoreMap(:,:,1) = sidsam(dcube, phaseLibrary(1).Reflectance);
scoreMap(:,:,2) = sidsam(dcube, phaseLibrary(2).Reflectance);
scoreMap(:,:,3) = sidsam(dcube, phaseLibrary(3).Reflectance);
scoreMap(:,:,4) = sidsam(dcube, phaseLibrary(4).Reflectance);
scoreMap(:,:,5) = sidsam(dcube, phaseLibrary(5).Reflectance);
scoreMap(:,:,6) = sidsam(dcube, phaseLibrary(6).Reflectance);

%% B4C
figure
imagesc(scoreMap(:,:,1))
colormap(jet);
colorbar
clim([0 700])

%% B6.5C
figure
imagesc(scoreMap(:,:,2))
colormap(jet);
colorbar
clim([0 1000])

%% Si-BC
figure
imagesc(scoreMap(:,:,3))
colormap(jet);
colorbar
clim([0 1000])

%% Si
figure
imagesc(scoreMap(:,:,4))
colormap(jet);
colorbar
clim([0 300])

%% SiB6
figure
imagesc(scoreMap(:,:,5))
colormap(jet);
colorbar
clim([0 550])

%% B
figure
imagesc(scoreMap(:,:,6))
colormap(jet);
colorbar
clim([0 1000])

% figure
% montage(scoreMap,'Size',[2 3],'BorderSize',10)
% title('Score Map Obtained for Each Pure Spectrum','FontSize',14)
% colormap(jet);
% colorbar
% clim([0 1000])