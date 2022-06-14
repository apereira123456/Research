%% Clear
close all; clear all; clc;

%% Setup Parameters
% Truncation Range (Wavenumber)
trunc = [150, 1500];

% Step size in microns
step = 1;

% Origin
origin_index = 961;
%401 404
% B4C after origin_index
% SiC before origin_index

% Interval for moving average
interval = 51;

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

ij = x_size*y_size;

%% Specify analysis region
% Prompt user to select csv file from Raman Map to Line program
if x_size ~= 1 && y_size ~= 1
    [file, path] = uigetfile({'*.csv'},'Select index file',path);
    index_name = fullfile(path, file);
    index = table2array(readtable(index_name));
    index = index';
% Prompt user to specify start and end indices for line scan
else
    prompt = "Enter start index:";
    a = input(prompt);

    prompt = "Enter end index:";
    b = input(prompt);

    if a < b
        index = a:b;
    else
        index = b:a;
    end
end

%% Convert raw data to plottable format
map_data = cell(x_size,y_size,2);

[~,trunc_start] = min(abs(data(1:s_size,3) - trunc(2)));
[~,trunc_end] = min(abs(data(1:s_size,3) - trunc(1)));

trunc_size = trunc_end - trunc_start + 1;

start_offset = trunc_start - 1;
end_offset = s_size - trunc_end;

for j = 1:y_size
    for i = 1:x_size
        a = (i-1)*s_size + (j-1)*x_size*s_size + 1;
        a = a + start_offset;

        b = i*s_size + (j-1)*x_size*s_size;
        b = b - end_offset;

        map_data(i,j,1) = mat2cell(data(b:-1:a,3),trunc_size,1);
        map_data(i,j,2) = mat2cell(normalize(data(b:-1:a,4),1,'range'),trunc_size,1);
    end
end

%% Graph Labels
title_text = '\textbf{Doped Boron Carbide - 1080 cm$^{-1}$ Peak Shift}';
subtitle_text = '(Adjacent to Undoped Boron Carbide)';
x_text = '\textbf{Distance from Interface ($\mu$m)}';
y_text = '\textbf{Raman Shift (cm$^{-1}$)}';

%% Peak Shift
peak_loc = zeros(1,ij);

for n = index
    i = mod(n, x_size);
    if i == 0
        i = x_size;
    end
    j = ceil(n/x_size);

    [M, I] = max(map_data{i,j,2}, [], 1);
    peak_loc(1,n) = map_data{i,j,1}(start_index + I - 1,1);
end

% Map Scans: Correct distance for lines on an anlge
if x_size ~= 1 && y_size ~=1
    x_coord = [mod(index(1), x_size), mod(index(end), x_size)];
    y_coord = [ceil(index(1) / x_size), ceil(index(end) / x_size)];
    
    correction = x_coord == 0;
    x_coord(correction) = x_size;

    line_length = step * sqrt((x_coord(end) - x_coord(1))^2 + (y_coord(end) - y_coord(1))^2);
    
    Distance = linspace(1,line_length,length(index));

% Line Scans: Adjust distance based on origin
else
    Distance = step * (index - origin_index);
end

figure()
plot(Distance, peak_loc(1,index), 'o')
title(title_text, 'interpreter', 'latex')
xlabel(x_text, 'interpreter', 'latex')
ylabel(y_text, 'interpreter', 'latex')
xlim([Distance(1) Distance(end)])

moving_avg = movmean(peak_loc(1,index),interval);

figure()
hold on
plot(Distance, moving_avg, 'b', 'LineWidth', 2)
title(title_text, 'interpreter', 'latex', 'FontSize', 18)
%subtitle(subtitle_text, 'interpreter', 'latex', 'FontSize', 14)
xlabel(x_text, 'interpreter', 'latex', 'FontSize', 14)
ylabel(y_text, 'interpreter', 'latex', 'FontSize', 14)
xlim([Distance(1) Distance(end)])

%% Standard Deviation
% moving_std = movstd(peak_loc(1,index),interval);
% 
% top = moving_avg + moving_std;
% bottom = moving_avg - moving_std;
% 
% x_std = [Distance, fliplr(Distance)];
% y_std = [top, fliplr(bottom)];
% 
% fill(x_std, y_std, [0.5 0.5 0.5], 'FaceAlpha', 0.3);