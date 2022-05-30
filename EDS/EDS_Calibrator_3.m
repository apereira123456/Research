% Kent Christian 4/28/22
% Select B, C, and Si maps from EDS saved as .csv files.
% Estimates the true At.% of these elements based on a standard calibration

close all; clear all; clc;

% Parameters that can be changed
% ======================================
% 2 scaling factor functions are used to adjust intensity for C and Si
% based on B intensity
% each of the 2 scaling factors has 2 fitting parameters (a and m)
% cb_scaling_factor = a*exp(m .* cb_intensity_ratio);
cb_a = 0.83;    cb_m = 0.75;
sb_a = 0.63;    sb_m = -0.45;

% set the pixel bounds of the region we want to look at
% xmin = 300; xmax = 475;
% ymin = 570; ymax = 670;
full_map = 1;
xmin = 408; xmax = 1024;
ymin = 1; ymax = 324;
% ======================================

%% ask user to locate the 3 .csv data files
[filename, path] = uigetfile('*','Select EDS data files (.cvs)','MultiSelect','on');
number_of_files = length(filename);

%% iterate over all files chosen
for ix = 1:number_of_files

    % create path to file
    mapFile = strcat(path,filename{ix});

    % read the data file
    data = dlmread(mapFile,','); % loads data from the file into a matrix

    % bin n by n pixels to reduce resolution of the image. Default is 4
    y_size = size(data,1);
    x_size = size(data,2);
    cb_scaling_factor = 4; % bins are n by n pixel squares where n=scaling_factor
    binned_y_size = floor(y_size/cb_scaling_factor);
    binned_x_size = floor(x_size/cb_scaling_factor);
    binned_data = imresize(data, [binned_y_size binned_x_size]);

    if full_map
        ymin = 1; ymax = binned_y_size;
        xmin = 1; xmax = binned_x_size; 
    end

    % plot the binned data map
    plotnumber = 2*ix - 1;
    subplot(number_of_files, 2, plotnumber)
    imshow(binned_data)
    colormap(summer)
    title(filename{ix} + " - binned map data")
  
    % draw a rectangle around the selected region in the binned map
    xdisplace = xmax-xmin;
    ydisplace = ymax-ymin;
    rectangle('Position',[xmin ymin xdisplace ydisplace],'EdgeColor','r','LineWidth',1)
    axis on

    % extract data from the region as its own matrix
    region = binned_data(ymin:ymax, xmin:xmax);

    % line average the intensities over the length of the region
    line_average(ix,:) = mean(region);
    full_length_line_avg(ix,:) = mean(binned_data);

    % plot the line averaged intensity of the map
    subplot(number_of_files, 2, plotnumber+1)
    full_length_x = 1:binned_x_size;
    plot(full_length_x, full_length_line_avg(ix,:), '.')
    title(strcat(filename{ix},' -',' line average intensity - log(x) fit'))
    xline(xmin); xline(xmax);

%     % do a log curve fit over the line averaged intensity
%     x = xmin:xmax;
%     logx = log(x);
%     [fit, ~] = polyfit(logx,line_average(ix,:),1);
%     m = fit(1);
%     b = fit(2);
%     lin_fit = m*logx+b;

%     % plot the curve fit over the line averaged intensity
%     subplot(number_of_files,2,plotnumber+1)
%     hold on
%     plot(x,lin_fit, 'LineWidth',1)
    
end

%% Apply the Scaling Factors
% save each element's line averaged intensity as its own vector
region_b_counts = line_average(1,:);
region_c_counts = line_average(2,:);
region_si_counts = line_average(3,:);

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

% calculate average atomic percents over the region
avg_at_percent_c = mean(at_percent_c);
avg_at_percent_si = mean(at_percent_si);
avg_at_percent_b = mean(at_percent_b);
b_to_c_ratio = avg_at_percent_b/avg_at_percent_c;

% print the average atomic percents
fprintf('\nAverage atomic %% over the region:\n')
fprintf('carbon %0.2f %%\n',avg_at_percent_c)
fprintf('silicon %0.2f %%\n',avg_at_percent_si)
fprintf('boron %0.2f %%\n',avg_at_percent_b)
fprintf('B to C ratio: %0.1f\n\n',b_to_c_ratio)

% plot the line averaged atomic percents over the length of the region
figure();
title('Atomic % over the EDS Region')
xlabel('Position (pixel)'), ylabel('Atomic Percent (%)')
hold on
plot(at_percent_c,'.')
plot(at_percent_si,'.')
plot(at_percent_b,'.')
legend('C','Si','B')
