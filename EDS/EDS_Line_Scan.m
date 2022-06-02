close all; clear all; clc;

% Setup parameters
full_map = 0;

%% Specify scaling factors depending on SEM
SEM = 'Gemini';

if strcmp(SEM,'Gemini')
    cb_a = 0.83;    cb_m = 0.75;
    sb_a = 0.63;    sb_m = -0.45;
elseif strcmp(SEM,'Sigma')
    cb_a = 0.83;    cb_m = 0.75;
    sb_a = 0.63;    sb_m = -0.45;
end

%% Prompt user to select the 3 data files
[filename, path] = uigetfile('*.csv','Select EDS data files','MultiSelect','on');
number_of_files = length(filename);

%% Read Data
map_size = size(readmatrix(strcat(path,filename{1})));
data = zeros(map_size(1),map_size(2),number_of_files);

for i = 1:number_of_files
    % create path to file
    mapFile = strcat(path,filename{i});

    % read the data file
    data(:,:,i) = readmatrix(mapFile); % loads data from the file into a matrix
end

[data_y_size, data_x_size, ~] = size(data);

%% Display EDS Map
figure()
hold on
imshow(data)
colormap(summer)

%% Specify analysis region
if full_map == 1
    xmin = 1; xmax = data_x_size;
    ymin = 1; ymax = data_y_size;

    w = xmax - xmin + 1;
    h = ymax - ymin + 1;

    rectangle('Position', [xmin ymin w h], 'EdgeColor', 'r', 'LineWidth', 2)
else
    [x_coord, y_coord] = ginput(2);
    
    x1 = x_coord(1); x2 = x_coord(2);
    y1 = y_coord(1); y2 = y_coord(2);
    
    % Determine corner and dimensions
    x = min(x1, x2);
    y = min(y1, y2);
    w = abs(x2 - x1);
    h = abs(y2 - y1);
    
    % Plot rectangle on top of image
    rectangle('Position', [x y w h], 'EdgeColor', 'r', 'LineWidth', 2)

    % Store values
    % Start is bottom left and end is top right
    xmin = floor(x); xmax = floor(x + w);
    ymin = floor(y); ymax = floor(y + h);
    w = xmax - xmin + 1;
end

[region_b_counts, region_c_counts, region_si_counts] = deal(zeros(1,w));

for i = 1:w
    region_b_counts(i) = mean(data(ymin:ymax, xmin + i - 1, 1), 'all');
    region_c_counts(i) = mean(data(ymin:ymax, xmin + i - 1, 2), 'all');
    region_si_counts(i) = mean(data(ymin:ymax, xmin + i - 1, 3), 'all');
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

% calculate average atomic percents over the region
avg_at_percent_c = mean(at_percent_c);
avg_at_percent_si = mean(at_percent_si);
avg_at_percent_b = mean(at_percent_b);
avg_b_to_c_ratio = mean(b_to_c_ratio);

% print the average atomic percents
fprintf('Average values over the region:\n')
% fprintf('B: %0.2f %%\n',avg_at_percent_b)
% fprintf('C: %0.2f %%\n',avg_at_percent_c)
fprintf('Si: %0.2f at.%%\n',avg_at_percent_si)
fprintf('B/C: %0.2f\n',avg_b_to_c_ratio)

% plot the line averaged atomic percents over the length of the region
figure()
plot(at_percent_si,'.')
title('Atomic % Si Over the EDS Region')
xlabel('Position (pixel)')
ylabel('Atomic Percent (%)')

figure()
plot(b_to_c_ratio,'.')
title('B/C Ratio Over the EDS Region')
xlabel('Position (pixel)')
ylabel('B/C Ratio')