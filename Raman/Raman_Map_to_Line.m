%% Clear
close all; clear all; clc;

%% Data Selection Prompt
[data_file, data_path] = uigetfile({'*.txt'}, 'Select map data file');
data_name = fullfile(data_path, data_file);
data = table2array(readtable(data_name));

%% Logic for Scans of Different Dimensions
% Determine data array size
[~,x,~] = unique(data(:,1),'rows');
[~,y,~] = unique(data(:,2),'rows');
[~,s,~] = unique(data(:,3),'rows');

[x_size, ~] = size(x);
[y_size, ~] = size(y);
[s_size, ~] = size(s);

%% Image selection prompt
[image_file, image_path] = uigetfile({'*.tif';'*.tiff';'*.png';'lpg';'*.jpeg'},'Select image of map region',data_path);
image_name = fullfile(image_path, image_file);
image = imread(image_name);

%% Coordinate selection prompt
figure()
hold on
imshow(image)
[x_coord, y_coord] = ginput(2);
line(x_coord, y_coord, 'Color', 'r', 'LineWidth', 2)

%%
% Line Intersections
% Determine image size and generate scale factor
[y_pixels, x_pixels, ~] = size(image);
scale = [x_pixels/(x_size), y_pixels/(y_size)];

% Scale x-values and determine min/max
x_pos = x_coord / scale(1);
x_min = min(x_pos, [], 1);
x_max = max(x_pos, [], 1);

% This step may not be necessary because the y-axis on the raman is flipped?
y_coord = y_pixels - y_coord;

% Scale y-values and determine min/max
y_pos = y_coord / scale(2);
y_min = min(y_pos, [], 1);
y_max = max(y_pos, [], 1);

% Create a grid from the min/max values
x = floor(x_min) : 1 : ceil(x_max);
y = floor(y_min) : 1 : ceil(y_max);

% Determine line equation from selected coordinates
m = (y_pos(2) - y_pos(1)) / (x_pos(2) - x_pos(1));
b = (y_pos(1) - m * x_pos(1));

hrz = @(y, m, b) [(y - b)./m; y];
vrt = @(x, m, b) [x; m.*x + b];

% Round up or down intercepts based on slope
if m < 0
    h_int = ceil(hrz(y, m, b)');
    v_int = floor(vrt(x, m, b)');
    hv_int = [h_int; v_int];
else
    h_int = ceil(hrz(y, m, b)');
    v_int = ceil(vrt(x, m, b)');
    hv_int = [h_int; v_int];
end

% Remove values that lie outside the grid boundary
ex_bdry = find(hv_int(:,1) < floor(x_min) | hv_int(:,1) < 1);
ex_bdry = [find(hv_int(:,1) > ceil(x_max) | hv_int(:,1) > x_size); ex_bdry];
ex_bdry = [find(hv_int(:,2) < y_min | hv_int(:,2) < 1 ); ex_bdry];
ex_bdry = [find(hv_int(:,2) > y_max | hv_int(:,2) > y_size); ex_bdry];
hv_int(ex_bdry,:) = [];

% Remove duplicate values
if m < 0
    unq = unique(hv_int,'rows');
    srtd = flip(sortrows(unq,[1 2],{'descend' 'ascend'}));
else
    srtd = unique(hv_int,'rows');
end

%% Write index values to csv file
index = srtd(:,1) + x_size * (srtd(:,2) - 1);
writematrix(index, fullfile(data_path,'Index.csv'));

% i = srtd(:,1);
% j = srtd(:,2);

a = [x_pos(1), x_pos(2)];
b = [y_pos(1), y_pos(2)];

%% Draw line and grid
figure()
plot(a,b)
xlim([1 x_size])
ylim([1 y_size])
set(gca,'xtick',1:1:x_size)
set(gca,'ytick',1:1:y_size)
grid on