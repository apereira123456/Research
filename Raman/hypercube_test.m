%% Clear
close all; clear all; clc;

%% Setup Parameters
% Specify excitation laser wavelength in nm
lambda_excitation = 532;

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
% Determine the indices associated with the truncation start and stop
trunc_range = [150, 1500];
[~,trunc_start] = min(abs(data(1:s_size,3) - trunc_range(2)));
[~,trunc_end] = min(abs(data(1:s_size,3) - trunc_range(1)));

trunc_size = trunc_end - trunc_start + 1;

% Adjust the truncation start and stop so the values can be used with indexing variables a and b
start_offset = trunc_start - 1;
end_offset = s_size - trunc_end;

a = (1-1)*s_size + (1-1)*x_size*s_size + 1;
a = a + start_offset;

b = 1*s_size + (1-1)*x_size*s_size;
b = b - end_offset;

wavelength = data(b:-1:a,3);

% Convert wavenumber (cm^-1) to wavelength (nm)
wavelength = 10000000 ./ wavelength;

% Initialize the cell that stores the Raman data
dcube = zeros(y_size,x_size,trunc_size);

for j = 1:y_size
    for i = 1:x_size
        a = (i-1)*s_size + (j-1)*x_size*s_size + 1;
        a = a + start_offset;

        b = i*s_size + (j-1)*x_size*s_size;
        b = b - end_offset;

        dcube(j,i,:) = normalize(data(b:-1:a,4),1,'range');
    end
end

hcube = hypercube(dcube(:,:,1:end-1), wavelength(1:end-1));
% hyperspectralViewer(hcube);

% newhcube = removeBands(hcube,'Wavelength',wlrange)
% newhcube = removeBands(hcube,'Wavelength', [-1000 150; 1500 3000]);

%% Standard Spectra
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
standard(:,1) = wavelength(1:end-1);

libData = struct('Name', {}, 'Type', {}, 'Class', {}, 'SubClass', {}, 'ParticleSize', {}, ...
                 'Genus', {}, 'Species', {}, 'SampleNo', {}, 'Owner', {}, 'WavelengthRange', {}, ...
                 'Origin', {}, 'CollectionDate', {}, 'Description', {}, 'Measurement', {}, ...
                 'FirstColumn', {}, 'SecondColumn', {}, 'WavelengthUnit', {}, 'DataUnit', {}, ...
                 'FirstXValue', {}, 'SecondXValue', {}, 'NumberOfXValues', {}, ...
                 'AdditionalInformation', {}, 'Wavelength', {}, 'Reflectance', {});

libData(1).Name = 'Boron Carbide'; libData(2).Name = 'Silicon'; libData(3).Name = 'Silicon Boride';
libData(1).Type = 'Ceramic'; libData(2).Type = 'Semiconductor'; libData(3).Type = 'Ceramic';
libData(1).Class= 'Starck'; libData(2).Class = 'Fisher'; libData(3).Class = 'Materion';
libData(1).SubClass = 'HD20'; libData(2).SubClass = '99.5%'; libData(3).SubClass = '98%';
libData(1).ParticleSize = '0.5um'; libData(2).ParticleSize = '-325 Mesh'; libData(3).ParticleSize = '-200 Mesh';
[libData(1).Genus, libData(2).Genus, libData(3).Genus] = deal('');
[libData(1).Species, libData(2).Species, libData(3).Species] = deal('');
libData(1).SampleNo = 'Reference B4C'; libData(2).SampleNo = 'Reference Si'; libData(3).SampleNo = 'Reference SiB6';
[libData(1).Owner, libData(2).Owner, libData(3).Owner] = deal('Dept. of Materials Science and Engineering, Rutgers University');
[libData(1).WavelengthRange, libData(2).WavelengthRange, libData(3).WavelengthRange] = deal('All');
[libData(1).Origin, libData(2).Origin, libData(3).Origin] = deal('Andrew Pereira');
[libData(1).CollectionDate, libData(2).CollectionDate, libData(3).CollectionDate] = deal('N/A');
[libData(1).Description, libData(2).Description, libData(3).Description] = deal('Reference spectrum for phase ID');
[libData(1).Measurement, libData(2).Measurement, libData(3).Measurement] = deal('Renishaw inVia Confocal Raman Microscope');
[libData(1).FirstColumn, libData(2).FirstColumn, libData(3).FirstColumn] = deal('X');
[libData(1).SecondColumn, libData(2).SecondColumn, libData(3).SecondColumn] = deal('Y');
[libData(1).WavelengthUnit, libData(2).WavelengthUnit, libData(3).WavelengthUnit] = deal('nanometer');
[libData(1).DataUnit, libData(2).DataUnit, libData(3).DataUnit] = deal('Intensity (counts)');
[libData(1).FirstXValue, libData(2).FirstXValue, libData(3).FirstXValue] = deal(66348.7455940286);
[libData(1).SecondXValue,libData(2).SecondXValue, libData(3).SecondXValue] = deal(6664.67941033257);
[libData(1).NumberOfXValues, libData(2).NumberOfXValues, libData(3).NumberOfXValues] = deal(486);
[libData(1).AdditionalInformation, libData(2).AdditionalInformation, libData(3).AdditionalInformation] = deal('none');
[libData(1).Wavelength, libData(2).Wavelength, libData(3).Wavelength] = deal(standard(:,1));
libData(1).Reflectance = standard(:,2); libData(2).Reflectance = standard(:,3); libData(3).Reflectance = standard(:,4);

%% Match Spectra
% figure
% ax(1) = subplot(4,1,1); 
% plot(wavelength, squeeze(dcube(1,1,:)))
% ylabel('Intensity')
% grid on
% title('Sample B4C Spectrum')
% 
% ax(2) = subplot(4,1,2); 
% plot(libData(1).Wavelength, libData(1).Reflectance)
% ylabel('Intensity') 
% grid on
% title('Reference B4C Spectrum')
% 
% ax(2) = subplot(4,1,3); 
% plot(libData(1).Wavelength, libData(2).Reflectance)
% ylabel('Intensity') 
% grid on
% title('Reference Si Spectrum')
% 
% ax(2) = subplot(4,1,4); 
% plot(libData(1).Wavelength, libData(3).Reflectance)
% ylabel('Intensity') 
% grid on
% title('Reference SiB6 Spectrum')
% 
% figure
% ax(1) = subplot(4,1,1); 
% plot(wavelength, squeeze(dcube(30,40,:)))
% ylabel('Intensity')
% grid on
% title('Sample B4C Spectrum')
% 
% ax(2) = subplot(4,1,2); 
% plot(libData(1).Wavelength, libData(1).Reflectance)
% ylabel('Intensity') 
% grid on
% title('Reference B4C Spectrum')
% 
% ax(2) = subplot(4,1,3); 
% plot(libData(1).Wavelength, libData(2).Reflectance)
% ylabel('Intensity') 
% grid on
% title('Reference Si Spectrum')
% 
% ax(2) = subplot(4,1,4); 
% plot(libData(1).Wavelength, libData(3).Reflectance)
% ylabel('Intensity') 
% grid on
% title('Reference SiB6 Spectrum')

scoreMap(:,:,1) = spectralMatch(libData(1),hcube);
scoreMap(:,:,2) = spectralMatch(libData(2),hcube);
scoreMap(:,:,3) = spectralMatch(libData(3),hcube);

figure
montage(scoreMap,'Size',[1 numel(libData)],'BorderSize',10)
title('Score Map Obtained for Each Pure Spectrum','FontSize',14)
colormap(jet);
colorbar
clim([0 5])

% scoreMap(:,:,1) = sidsam(dcube(:,:,1:end-1), libData(1).Reflectance);
% scoreMap(:,:,2) = sidsam(dcube(:,:,1:end-1), libData(2).Reflectance);
% scoreMap(:,:,3) = sidsam(dcube(:,:,1:end-1), libData(3).Reflectance);
% 
% figure
% montage(scoreMap,'Size',[1 numel(libData)],'BorderSize',10)
% title('Score Map Obtained for Each Pure Spectrum','FontSize',14)
% colormap(jet);
% colorbar
% clim([0 100])