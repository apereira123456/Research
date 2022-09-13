%% Standard Spectra
% Read standard spectra data into array
standard = readtable('Standard Raman Spectra.txt');

labels = standard.Properties.VariableNames;

labels(2:end)

% standard = flipud(standard);