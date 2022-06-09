close all; clear all; clc;

% Setup parameters
full_map = 0;

% Atomic masses: B, C, Si (g/mol)
atomic_mass = [10.811, 12.011, 28.0855];

%% Specify scaling factors depending on SEM
SEM = 'Gemini';
% SEM = 'Sigma';

if strcmp(SEM,'Gemini')
    m = [0.75, -0.45];
    b = [0.83, 0.63];
elseif strcmp(SEM,'Sigma')
    m = [0.75, -0.45];
    b = [0.83, 0.63];
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

%% Specify analysis region
if full_map == 1
    x = 1; y = 1;

    xmin = 1; xmax = data_x_size;
    ymin = 1; ymax = data_y_size;

    w = xmax - xmin;
    h = ymax - ymin;

    rectangle('Position', [x y w h], 'EdgeColor', 'r', 'LineWidth', 2)

    delta = w + 1;
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
    % Start is top left and end is bottom right (image origin is top left)
    xmin = ceil(x); xmax = floor(x + w);
    ymin = ceil(y); ymax = floor(y + h);
    delta = xmax - xmin + 1;
end

counts = zeros(3,delta);

for i = 1:delta
    counts(1,i) = mean(data(ymin:ymax, xmin + i - 1, 1));
    counts(2,i) = mean(data(ymin:ymax, xmin + i - 1, 2));
    counts(3,i) = mean(data(ymin:ymax, xmin + i - 1, 3));
end

%% Apply the Scaling Factors
for i = 1:2
    % compute elemental intensity ratios to be used in scaling factor calcs
    intensity_ratio(i,:) = counts(i+1,:) ./ counts(1,:);
    % calculate scaling factors for each point along line averaged intensity
    scaling_factor(i,:) = b(i) * exp(m(i) .* intensity_ratio(i,:));
    % apply scaling factors to get corrected relative intensity of each element
    rel_intensity(i+1,:) = scaling_factor(i,:) .* intensity_ratio(i,:);
end

rel_intensity(1,:) = ones(1,size(rel_intensity,2));

total_intensity = rel_intensity(1,:) + rel_intensity(2,:) + rel_intensity(3,:);
for i = 1:3
    % convert corrected intensities to weight fractions
    wt_fraction(i,:) = rel_intensity(i,:) ./ total_intensity;
    % convert weight fractions to relative mol fraction
    mol(i,:) = wt_fraction(i,:) / atomic_mass(i);
end

% convert relative mol fraction to atomic percents
total_mols = mol(1,:) + mol(2,:) + mol(3,:);
for i = 1:3
    at_percent(i,:) = 100 * mol(i,:) ./ total_mols;
end

for i = 1:3
    % calculate average atomic percents over the region
    avg_at_percent(i) = mean(at_percent(i,:),2);
end

avg_b_to_c_ratio = avg_at_percent(1) / avg_at_percent(2);

b_to_c_ratio = at_percent(1,:) ./ at_percent(2,:);
isinf(b_to_c_ratio)
b_to_c_ratio(isinf(b_to_c_ratio)) = NaN;
mean(b_to_c_ratio,'omitnan')

% print the average atomic percents
fprintf('Average values over the region:\n')
fprintf('B: %0.2f at.%%\n',avg_at_percent(1))
fprintf('C: %0.2f at.%%\n',avg_at_percent(2))
fprintf('Si: %0.2f at.%%\n',avg_at_percent(3))
fprintf('B/C: %0.2f\n',avg_b_to_c_ratio)

%% Individual Plots
interval = xmin:xmax;

for i = 1:3
    subplot(number_of_files, 2, 2*i-1)
    imshow(data(:,:,i))
    colormap(summer)
    title(filename{i})
    
    rectangle('Position', [x y w h], 'EdgeColor', 'r', 'LineWidth', 1)
    
    subplot(number_of_files, 2, 2*i)
    plot(interval, counts(i,:), '.')
    title(strcat(filename{i},' -',' Line Averaged Intensity '))
    xline(xmin); xline(xmax);
end

%% Plot the line averaged atomic percents over the length of the region
figure()
hold on
for i = 1:3
    plot(at_percent(i,:),'.')
end

title('Atomic % Over the EDS Region')
xlabel('Position (pixel)')
ylabel('Atomic Percent (%)')
legend('B','C','Si')

figure()
hold on

yyaxis left
plot(at_percent(3,:),'.')
ylabel('Atomic % Si')

yyaxis right
plot(b_to_c_ratio,'.')
ylabel('B/C Ratio')

title('Atomic % Si and B/C Ratio Over the EDS Region')
xlabel('Position (pixel)')
legend('Si','B/C')