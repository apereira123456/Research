%%
% close all; 
clear all; clc

% Setup Parameters
m = [0.75, -0.45];                          % Scaling factors for C/B and Si/B
b = [0.83, 0.63];                           % Scaling factors for C/B and Si/B
atomic_mass = [10.811, 12.011, 28.0855];    % Atomic mass of B, C and Si

%% Prompt user to select 3 data files
[filename, path] = uigetfile('*.csv', 'Select EDS data files', 'MultiSelect', 'on');

% Determine number of files selected
number_of_files = length(filename);

%% Read Data
[y_size, x_size] = size(readmatrix(strcat(path,filename{1})));
ij = x_size * y_size;

data = zeros(y_size,x_size,number_of_files);

for i = 1:number_of_files
    % Create path to file
    mapFile = strcat(path,filename{i});

    % Read the data file
    data(:,:,i) = readmatrix(mapFile);
end

% Find the average of an nxn box of pixels surrounding a central pixel and replace it with the average 
n = 49;
% n = 75;
% n = 99;
c = ceil(n/2);
kernel = ones(n);
kernel(c,c) = 0;
for i = 1:3
    data(:,:,i) = conv2(data(:,:,i), kernel, 'same');
end

%% Compute average intensity for each grain
[counts, intensity_ratio, scaling_factor, rel_intensity, wt_fraction, mol, at_percent] = deal(ones(ij,3));

for index = 1:ij
    i = mod(index, x_size);
    if i == 0
        i = x_size;
    end
    j = ceil(index / x_size);

    counts(index,1) = data(j,i,1);
    counts(index,2) = data(j,i,2);
    counts(index,3) = data(j,i,3);
end

%% Apply the Scaling Factors
for i = 2:3
    % Compute elemental intensity ratios to be used in scaling factor calcs
    intensity_ratio(:,i) = counts(:,i) ./ counts(:,1);
    % Calculate scaling factors for each point along line averaged intensity
    scaling_factor(:,i) = b(i-1) * exp(m(i-1) .* intensity_ratio(:,i));
    % Apply scaling factors to get corrected relative intensity of each element
    rel_intensity(:,i) = scaling_factor(:,i) .* intensity_ratio(:,i);
end

% Compute the total intensity
total_intensity = rel_intensity(:,1) + rel_intensity(:,2) + rel_intensity(:,3);

for i = 1:3
    % Convert corrected intensities to weight fractions
    wt_fraction(:,i) = rel_intensity(:,i) ./ total_intensity;
    % Convert weight fractions to relative mol fraction
    mol(:,i) = wt_fraction(:,i) / atomic_mass(i);
end

% Compute the total number of moles
total_mols = mol(:,1) + mol(:,2) + mol(:,3);

for i = 1:3
    % Convert relative mol fraction to atomic percents
    at_percent(:,i) = 100 * mol(:,i) ./ total_mols;
end

% % Compute the B to C ratio
% b_to_c_ratio = at_percent(:,1) ./ at_percent(:,2);

[boron, carbon, silicon] = deal(NaN(y_size,x_size));

for index = 1:ij
    i = mod(index, x_size);
    if i == 0
        i = x_size;
    end
    j = ceil(index / x_size);

    boron(j,i) = at_percent(index,1);
    carbon(j,i) = at_percent(index,2);
    silicon(j,i) = at_percent(index,3);
end

b_to_c_ratio = boron ./ carbon;

%%
figure
silicon_map = imagesc(silicon);
set(silicon_map,'AlphaData',~isnan(silicon))
colormap jet
colorbar
clim([1.7 3.2]);

%%
figure
bc_map = imagesc(b_to_c_ratio);
set(bc_map,'AlphaData',~isnan(b_to_c_ratio))
colormap jet
colorbar
clim([0 3.5]);