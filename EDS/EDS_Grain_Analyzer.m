close all; clear all; clc;

% Parameters that can be changed
% ======================================
% 2 scaling factor functions are used to adjust intensity for C and Si
% based on B intensity
% each of the 2 scaling factors has 2 fitting parameters (a and m)
% cb_scaling_factor = a*exp(m .* cb_intensity_ratio);
cb_a = 0.83;    cb_m = 0.75;
sb_a = 0.63;    sb_m = -0.45;

%% Ask user to locate the 3 .csv data files
[filename, path] = uigetfile('*','Select EDS data files (.cvs)','MultiSelect','on');
number_of_files = length(filename);


%% Plot Data
for i = 1:number_of_files
    % create path to file
    mapFile = strcat(path,filename{i});

    % read the data file
    data(:,:,i) = dlmread(mapFile,','); % loads data from the file into a matrix
    
%     figure()
%     imshow(data(:,:,i))
end

[data_y_size, data_x_size, ~] = size(data);

%% Grain selection prompt
prompt = "Enter number of grains:";
num_grains = input(prompt);

[x_start, x_end, y_start, y_end] = deal(zeros(num_grains, 1));

figure()
hold on
imshow(data)

for i = 1:num_grains
    % Choose corners of rectangle
    [x_coord, y_coord] = ginput(2);

    x1 = x_coord(1); x2 = x_coord(2);
    y1 = y_coord(1); y2 = y_coord(2);

    x = min(x1, x2);
    y = min(y1, y2);
    w = abs(x2 - x1);
    h = abs(y2 - y1);
    
    % Plot rectangle on top of image
    rectangle('Position', [x y w h], 'EdgeColor', 'r', 'LineWidth', 2)

    % Store values
    % Start is bottom left and end is top right
    x_start(i) = floor(x); x_end(i) = floor(x + w);
    y_start(i) = floor(y); y_end(i) = floor(y + h);
end

%% Compute average intensity for each grain
for i = 1:num_grains
    region_b_counts(i) = mean(data(y_start(i):y_end(i), x_start(i):x_end(i), 1), 'all');
    region_c_counts(i) = mean(data(y_start(i):y_end(i), x_start(i):x_end(i), 2), 'all');
    region_si_counts(i) = mean(data(y_start(i):y_end(i), x_start(i):x_end(i), 3), 'all');
end

%% Apply the Scaling Factors

% compute elemental intensity ratios to be used in scaling factor calcs
cb_intensity_ratio = region_c_counts ./ region_b_counts;
sib_intensity_ratio = region_si_counts ./ region_b_counts;

% calculate scaling factors for each point along line averaged intensity
cb_scaling_factor = 0.83*exp(0.75 .* cb_intensity_ratio);
sb_scaling_factor = 0.63*exp(-0.45 .* sib_intensity_ratio);

% apply scaling factors to get corrected relative intensity of each element
rel_intensity_c = cb_scaling_factor .* cb_intensity_ratio;
rel_intensity_si = sb_scaling_factor .* sib_intensity_ratio;
rel_intensity_b = ones(1,size(region_b_counts,2)); % B is 1 as the reference

% convert corrected intensities to weight fractions
total_intensity = rel_intensity_c + rel_intensity_si + rel_intensity_b;
wt_fraction_c = rel_intensity_c ./ total_intensity;
wt_fraction_si = rel_intensity_si ./ total_intensity;
wt_fraction_b = rel_intensity_b ./ total_intensity;

% atomic masses
c_mass = 12.011; % grams/mol
si_mass = 28.0855; % grams/mol
b_mass = 10.811; % grams/mol

% convert weight fractions to relative mol fraction
mol_c = wt_fraction_c / c_mass;
mol_si = wt_fraction_si / si_mass;
mol_b = wt_fraction_b / b_mass;

% convert relative mol fraction to atomic percents
total_mols = mol_c + mol_si + mol_b;
at_percent_c = 100 * mol_c ./ total_mols;
at_percent_si = 100 * mol_si ./ total_mols;
at_percent_b = 100 * mol_b ./ total_mols;
b_to_c_ratio = at_percent_b ./ at_percent_c;

%% print the average atomic percents
for i = 1:num_grains
    fprintf('Grain %i: \n', i)
    fprintf('Silicon: %0.2f at. %% \n', at_percent_si(i))
    fprintf('B to C ratio: %0.2f \n\n', b_to_c_ratio(i))
end