%% Clear
close all; clear all; clc;

%% Setup Parameters
% Specify excitation laser wavelength in nm
lambda_excitation = 532;
wave_start = 535.6;
wave_end = 595;

%% Standard Spectra
% Read standard spectra data into array
standard = table2array(readtable('Standard Raman Spectra.txt','ReadVariableNames',false));

% Sort the data in ascending order by Raman shift wavenumber
standard = flipud(standard);

% Convert Raman shift to Raman wavelength
standard(:,1) = ( (lambda_excitation)^(-1) - (standard(:,1) ./ 10^7) ).^(-1);

[~, start_index] = min(abs(standard(:,1) - wave_start));
[~, end_index] = min(abs(standard(:,1) - wave_end));

standard = standard(start_index:end_index, :);

% Parameters for library creation
[NumberOfXValues, NumberOfStandards] = size(standard);

FirstXValue = standard(1,1);
LastXValue = standard(end,1);

standard(:,2:NumberOfStandards) = normalize(standard(:,2:NumberOfStandards),1,'range');

%% Plot Standard Spectra
NumberOfStandards = NumberOfStandards - 1;

figure
for i = 1:NumberOfStandards
    ax(i) = subplot(NumberOfStandards,1,i); 
    plot(standard(:,1), standard(:,i + 1))
    grid on
end

%% Library Creation
phaseLibrary = struct('Name', {}, 'Type', {}, 'Class', {}, 'SubClass', {}, ...
                      'ParticleSize', {}, 'Genus', {}, 'Species', {}, ...
                      'SampleNo', {}, 'Owner', {}, 'WavelengthRange', {}, ...
                      'Origin', {}, 'CollectionDate', {}, 'Description', {}, ...
                      'Measurement', {}, 'FirstColumn', {}, 'SecondColumn', {}, ...
                      'WavelengthUnit', {}, 'DataUnit', {}, 'FirstXValue', {}, ...
                      'LastXValue', {}, 'NumberOfXValues', {}, ...
                      'AdditionalInformation', {}, 'Wavelength', {}, ...
                      'Reflectance', {});

phaseLibrary(1).Name = "Boron Carbide";
phaseLibrary(1).Type = "Ceramic";
phaseLibrary(1).Class = "Starck";
phaseLibrary(1).SubClass = "HD20";
phaseLibrary(1).ParticleSize = "0.5um";
phaseLibrary(1).Genus = "";
phaseLibrary(1).Species = "";
phaseLibrary(1).SampleNo = "Reference B4C";
phaseLibrary(1).Owner = "Dept. of Materials Science and Engineering, Rutgers University";
phaseLibrary(1).WavelengthRange = "All";
phaseLibrary(1).Origin = "Andrew Pereira";
phaseLibrary(1).CollectionDate = "N/A";
phaseLibrary(1).Description = "Reference spectrum for phase ID";
phaseLibrary(1).Measurement = "Renishaw inVia Confocal Raman Microscope";
phaseLibrary(1).FirstColumn = "X";
phaseLibrary(1).SecondColumn = "Y";
phaseLibrary(1).WavelengthUnit = "nanometer";
phaseLibrary(1).DataUnit = "Intensity (counts)";
phaseLibrary(1).FirstXValue = FirstXValue;
phaseLibrary(1).LastXValue = LastXValue;
phaseLibrary(1).NumberOfXValues = NumberOfXValues;
phaseLibrary(1).AdditionalInformation = "45W 532nm laser - 1200l/mm grating - 1500cm^-1 center - 100X objective - 1s exposure - 50% laser power - 50 accumulations";
phaseLibrary(1).Wavelength = standard(:,1);
phaseLibrary(1).Reflectance = standard(:,2);

phaseLibrary(2).Name = "Silicon";
phaseLibrary(2).Type = "Metalloid";
phaseLibrary(2).Class = "Fisher";
phaseLibrary(2).SubClass = "99.5%";
phaseLibrary(2).ParticleSize = "-325 Mesh";
phaseLibrary(2).Genus = "";
phaseLibrary(2).Species = "";
phaseLibrary(2).SampleNo = "Reference Si";
phaseLibrary(2).Owner = "Dept. of Materials Science and Engineering, Rutgers University";
phaseLibrary(2).WavelengthRange = "All";
phaseLibrary(2).Origin = "Andrew Pereira";
phaseLibrary(2).CollectionDate = "N/A";
phaseLibrary(2).Description = "Reference spectrum for phase ID";
phaseLibrary(2).Measurement = "Renishaw inVia Confocal Raman Microscope";
phaseLibrary(2).FirstColumn = "X";
phaseLibrary(2).SecondColumn = "Y";
phaseLibrary(2).WavelengthUnit = "nanometer";
phaseLibrary(2).DataUnit = "Intensity (counts)";
phaseLibrary(2).FirstXValue = FirstXValue;
phaseLibrary(2).LastXValue = LastXValue;
phaseLibrary(2).NumberOfXValues = NumberOfXValues;
phaseLibrary(2).AdditionalInformation = "45W 532nm laser - 1200l/mm grating - 1500cm^-1 center - 100X objective - 1s exposure - 50% laser power - 50 accumulations";
phaseLibrary(2).Wavelength = standard(:,1);
phaseLibrary(2).Reflectance = standard(:,3);

phaseLibrary(3).Name = "Silicon Boride";
phaseLibrary(3).Type = "Ceramic";
phaseLibrary(3).Class = "Materion";
phaseLibrary(3).SubClass = "98%";
phaseLibrary(3).ParticleSize = "-200 Mesh";
phaseLibrary(3).Genus = "";
phaseLibrary(3).Species = "";
phaseLibrary(3).SampleNo = "Reference SiB6";
phaseLibrary(3).Owner = "Dept. of Materials Science and Engineering, Rutgers University";
phaseLibrary(3).WavelengthRange = "All";
phaseLibrary(3).Origin = "Andrew Pereira";
phaseLibrary(3).CollectionDate = "N/A";
phaseLibrary(3).Description = "Reference spectrum for phase ID";
phaseLibrary(3).Measurement = "Renishaw inVia Confocal Raman Microscope";
phaseLibrary(3).FirstColumn = "X";
phaseLibrary(3).SecondColumn = "Y";
phaseLibrary(3).WavelengthUnit = "nanometer";
phaseLibrary(3).DataUnit = "Intensity (counts)";
phaseLibrary(3).FirstXValue = FirstXValue;
phaseLibrary(3).LastXValue = LastXValue;
phaseLibrary(3).NumberOfXValues = NumberOfXValues;
phaseLibrary(3).AdditionalInformation = "45W 532nm laser - 1200l/mm grating - 1500cm^-1 center - 100X objective - 1s exposure - 50% laser power - 50 accumulations";
phaseLibrary(3).Wavelength = standard(:,1);
phaseLibrary(3).Reflectance = standard(:,4);

phaseLibrary(4).Name = "Silicon Carbide";
phaseLibrary(4).Type = "Ceramic";
phaseLibrary(4).Class = "Starck";
phaseLibrary(4).SubClass = "UF-25";
phaseLibrary(4).ParticleSize = "0.45um";
phaseLibrary(4).Genus = "";
phaseLibrary(4).Species = "";
phaseLibrary(4).SampleNo = "Reference SiC";
phaseLibrary(4).Owner = "Dept. of Materials Science and Engineering, Rutgers University";
phaseLibrary(4).WavelengthRange = "All";
phaseLibrary(4).Origin = "Andrew Pereira";
phaseLibrary(4).CollectionDate = "N/A";
phaseLibrary(4).Description = "Reference spectrum for phase ID";
phaseLibrary(4).Measurement = "Renishaw inVia Confocal Raman Microscope";
phaseLibrary(4).FirstColumn = "X";
phaseLibrary(4).SecondColumn = "Y";
phaseLibrary(4).WavelengthUnit = "nanometer";
phaseLibrary(4).DataUnit = "Intensity (counts)";
phaseLibrary(4).FirstXValue = FirstXValue;
phaseLibrary(4).LastXValue = LastXValue;
phaseLibrary(4).NumberOfXValues = NumberOfXValues;
phaseLibrary(4).AdditionalInformation = "45W 532nm laser - 1200l/mm grating - 1500cm^-1 center - 100X objective - 1s exposure - 50% laser power - 50 accumulations";
phaseLibrary(4).Wavelength = standard(:,1);
phaseLibrary(4).Reflectance = standard(:,5);

% phaseLibrary(2).Name = "Boron Doped Boron Carbide";
% phaseLibrary(2).Type = "Ceramic";
% phaseLibrary(2).Class = "N/A";
% phaseLibrary(2).SubClass = "N/A";
% phaseLibrary(2).ParticleSize = "N/A";
% phaseLibrary(2).Genus = "";
% phaseLibrary(2).Species = "";
% phaseLibrary(2).SampleNo = "Reference B6.5C";
% phaseLibrary(2).Owner = "Dept. of Materials Science and Engineering, Rutgers University";
% phaseLibrary(2).WavelengthRange = "All";
% phaseLibrary(2).Origin = "Andrew Pereira";
% phaseLibrary(2).CollectionDate = "N/A";
% phaseLibrary(2).Description = "Reference spectrum for phase ID";
% phaseLibrary(2).Measurement = "Renishaw inVia Confocal Raman Microscope";
% phaseLibrary(2).FirstColumn = "X";
% phaseLibrary(2).SecondColumn = "Y";
% phaseLibrary(2).WavelengthUnit = "nanometer";
% phaseLibrary(2).DataUnit = "Intensity (counts)";
% phaseLibrary(2).FirstXValue = FirstXValue;
% phaseLibrary(2).LastXValue = LastXValue;
% phaseLibrary(2).NumberOfXValues = NumberOfXValues;
% phaseLibrary(2).AdditionalInformation = "none";
% phaseLibrary(2).Wavelength = standard(:,1);
% phaseLibrary(2).Reflectance = standard(:,6);
% 
% phaseLibrary(3).Name = "Silicon Doped Boron Carbide";
% phaseLibrary(3).Type = "Ceramic";
% phaseLibrary(3).Class = "N/A";
% phaseLibrary(3).SubClass = "N/A";
% phaseLibrary(3).ParticleSize = "N/A";
% phaseLibrary(3).Genus = "";
% phaseLibrary(3).Species = "";
% phaseLibrary(3).SampleNo = "Reference Si-BC";
% phaseLibrary(3).Owner = "Dept. of Materials Science and Engineering, Rutgers University";
% phaseLibrary(3).WavelengthRange = "All";
% phaseLibrary(3).Origin = "Andrew Pereira";
% phaseLibrary(3).CollectionDate = "N/A";
% phaseLibrary(3).Description = "Reference spectrum for phase ID";
% phaseLibrary(3).Measurement = "Renishaw inVia Confocal Raman Microscope";
% phaseLibrary(3).FirstColumn = "X";
% phaseLibrary(3).SecondColumn = "Y";
% phaseLibrary(3).WavelengthUnit = "nanometer";
% phaseLibrary(3).DataUnit = "Intensity (counts)";
% phaseLibrary(3).FirstXValue = FirstXValue;
% phaseLibrary(3).LastXValue = LastXValue;
% phaseLibrary(3).NumberOfXValues = NumberOfXValues;
% phaseLibrary(3).AdditionalInformation = "none";
% phaseLibrary(3).Wavelength = standard(:,1);
% phaseLibrary(3).Reflectance = standard(:,7);
% 
% phaseLibrary(6).Name = "Beta Boron";
% phaseLibrary(6).Type = "Metalloid";
% phaseLibrary(6).Class = "Starck";
% phaseLibrary(6).SubClass = "I";
% phaseLibrary(6).ParticleSize = "1.0-2.0 um";
% phaseLibrary(6).Genus = "";
% phaseLibrary(6).Species = "";
% phaseLibrary(6).SampleNo = "Reference B";
% phaseLibrary(6).Owner = "Dept. of Materials Science and Engineering, Rutgers University";
% phaseLibrary(6).WavelengthRange = "All";
% phaseLibrary(6).Origin = "Andrew Pereira";
% phaseLibrary(6).CollectionDate = "N/A";
% phaseLibrary(6).Description = "Reference spectrum for phase ID";
% phaseLibrary(6).Measurement = "Renishaw inVia Confocal Raman Microscope";
% phaseLibrary(6).FirstColumn = "X";
% phaseLibrary(6).SecondColumn = "Y";
% phaseLibrary(6).WavelengthUnit = "nanometer";
% phaseLibrary(6).DataUnit = "Intensity (counts)";
% phaseLibrary(6).FirstXValue = FirstXValue;
% phaseLibrary(6).LastXValue = LastXValue;
% phaseLibrary(6).NumberOfXValues = NumberOfXValues;
% phaseLibrary(6).AdditionalInformation = "none";
% phaseLibrary(6).Wavelength = standard(:,1);
% phaseLibrary(6).Reflectance = standard(:,8);

%% Save File 
save('phaseLibrary.mat','phaseLibrary')