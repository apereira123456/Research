%% Clear
close all; clear all; clc;

%% Setup Parameters
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
sorted_data = cell(x_size,y_size,2);
sorted_data(:,:,1) = mat2cell(data(s_size:-1:1,3),s_size,1);
% sorted_data(:,:,1) = data(s_size:-1:1,3); 

for j = 1:y_size
    for i = 1:x_size
        a = (i-1)*s_size + (j-1)*x_size*s_size + 1;
        b = i*s_size + (j-1)*x_size*s_size;

        sorted_data(i,j,2) = mat2cell(normalize(data(b:-1:a,4),1,'range'),s_size,1);
    end
end

Wavenumber = sorted_data{1,1,1};

%% Display Range Prompt
check = 1;

start_wavenumber = Wavenumber(1,1);
end_wavenumber = Wavenumber(end,1);

fprintf('Start Wavenumber: %.2f \n',start_wavenumber)
fprintf('End Wavenumber: %.2f \n',end_wavenumber)
    
while check
    prompt = "Enter wavenumber start value:";
    a = input(prompt);
    
    % Use the full range if no input is given
    if isempty(a)
        start_index = 1;
        end_index = s_size;
        disp_range = start_index:end_index;
        check = 0;
    
    % Otherwise find the index of the closest wavenumber
    else
        [~,start_index] = min(abs(Wavenumber - a));
        prompt = "Enter wavenumber end value:";
        a = input(prompt);
        [~,end_index] = min(abs(Wavenumber - a));
    
        disp_range = start_index:end_index;
        
        % Check for valid input
        if start_index > end_index
            disp('Invalid input.')

        elseif length(disp_range) > 1
            check = 0;

        else
            disp('Invalid input.')
        end
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

    [M, I] = max(sorted_data{i,j,2}(disp_range,1), [], 1);
    peak_loc(1,n) = sorted_data{i,j,1}(start_index + I - 1,1);
end

% Map Scans: Correct distance for lines on an anlge
if x_size ~= 1 && y_size ~=1
    x_coord = [mod(index(1), x_size), mod(index(end), x_size)];
    y_coord = [ceil(index(1)/x_size), ceil(index(end)/x_size)];
    
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