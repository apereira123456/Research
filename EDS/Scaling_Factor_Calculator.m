% Andrew Pereira 01/09/2023

close all; clear all; clc;

%% Locate the 3 .csv data files generated on the EDS to be calibrated
[filename, path] = uigetfile('*.csv','Select EDS data files (.cvs)','MultiSelect','on');
number_of_files = length(filename);

%% get the EDS data and refine it
% Parameters that can be changed
% ======================================
% set the pixel bounds of the region we want to look at
scale_factor = 1 / 4;
pixels_per_micron = scale_factor / 0.08390;
xmin = ceil(35 * pixels_per_micron); xmax = floor(xmin + 285 * pixels_per_micron);

SiB_intensity = [sum([81281.1,78790.3,79918.1,82409.5,81281.1]);
                 sum([31769.7,32062.2,32866.7,32755.3,31769.7]);
                 sum([88687.0,85556.6,83804.1,79707.1,88687.0])];

B4C_intensity = [sum([171121.4,184183.0,181401.9,159297.8,146274.2]);
                 sum([79732.9,79574.6,75086.0,73214.1,77753.7])];
% ymin = 1; y_max
% ======================================

% iterate over all files chosen
figure
for i = 1:number_of_files

    % create path to file
    mapFile = strcat(path,filename{i});

    % read the data file
    data = dlmread(mapFile,','); % loads data from the file into a matrix
    data = imresize(data, scale_factor, 'bicubic');

    y_size = size(data,1);
    x_size = size(data,2);

    % plot the data map
    subplot(2, 2, i)
    imshow(data)
    colormap(gray)
    title(filename{i} + " - map data")
  
    ymax = y_size;
    ymin = ceil(ymax - 215 * pixels_per_micron);

    % draw a rectangle around the selected region in the binned map
    xdisplace = xmax - xmin; % width of region
    ydisplace = ymax - ymin; % height of region
    rectangle('Position',[xmin ymin xdisplace ydisplace],'EdgeColor','r','LineWidth',1)
    axis on

    % extract data from the region as its own matrix
    region = data(ymin:ymax, xmin:xmax);

    % line average the intensities over the length of the region
    line_integrate(i,:) = sum(region,1);    
end

% Extra code if oxygen is plotted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
line_integrate(3,:) = [];
number_of_files = number_of_files - 1;
filename{1,3} = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Compute elemental intensity ratios to be used in scaling factor calcs
% Chris recommended integrating intensities, but since we take the ratio it
% should be the same either way
intensity_ratios = line_integrate ./ line_integrate(1,:); % divide each intensity by B's intensity

% Define regions to average intensities over. The first row is the start position
% and the second row is the end position in microns.
tem_regions = [0.01, 15, 85, 125, 165, 195, 240, 275;
               10, 25, 95, 135, 175, 205, 250, 285];

% Average the intensity ratios within defined regions
for i = 1:length(tem_regions)
    region_start = ceil(pixels_per_micron*tem_regions(1,i));
    region_end = floor(pixels_per_micron*tem_regions(2,i));
    average_tem_regions(:,i) = mean(intensity_ratios(:,region_start:region_end),2);
end

average_tem_regions(:,1) = SiB_intensity ./ SiB_intensity(1);
average_tem_regions(2,end) = B4C_intensity(2) ./ B4C_intensity(1);

% Standard data provided by Chris Marvel
at_percent = [86.1, 85.2, 84.9, 85.0, 84.7, 84.9, 84.8, 81.3;
              12.2, 12.8, 13.3, 13.1, 13.7, 13.8, 14.5, 18.6;
              1.8, 2.0, 1.8, 1.8, 1.6, 1.3, 0.7, 0.1];

% B, C, Si atomic masses
atomic_masses = [10.811 12.011 28.0855];

% Convert atomic percent to weight
for i = 1:number_of_files
    wt(i,:) = at_percent(i,:) * atomic_masses(i);
end

% Calculate total weight
total_wt = sum(wt,1);

% Calculate weight fraction
for i = 1:number_of_files
    wt_fraction(i,:) = wt(i,:) ./ total_wt;
end

% Convert weight fraction to weight ratio with respect to boron
wt_ratios = wt_fraction ./ wt_fraction(1,:);

% Calculate ZAFs at the TEM regions across standard
ZAF = wt_ratios ./ average_tem_regions; % divide by since average_tem_regions is I/I_B and we want I_B/I

name = {'Boron', 'Carbon', 'Silicon'};
% First guess for exponential fit based on previous fitting parameters
first_guess = [1, 1; 0.83, 0.75; 0.63, -0.45];

%% Plot ZAFs vs Intensity Ratios
figure
for i = 2:number_of_files
    subplot(1, number_of_files-1, i-1);

    hold on
    title(name(i))
    xlabel('Intensity Ratio') 
    ylabel('Scaling Factor')
    
    % Store intensity ratios as x (first point excluded in pure SiB6
    % because it is NaN)
    x = average_tem_regions(i,:);

    % Define functional form of equation we want to fit
    f = @(b,x) b(1).*exp(b(2).*x);

    % Optimize the fitting parameters, starting at first guess from
    % original fitting parameters
    B = fminsearch(@(b) norm(ZAF(i,:) - f(b,x)), first_guess(i,:));
    
    % Plot data and fit
    delta = max(x) - min(x);
    test = min(x) - delta:0.01:max(x) + delta;
    scatter(average_tem_regions(i,:), ZAF(i,:))
    plot(test, f(B,test))
    
    % Display equation
    legend(sprintf('f(x) = %.3f\\cdote^{%.3f\\cdotx}', B))
end