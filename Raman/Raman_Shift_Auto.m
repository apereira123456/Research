%% Clear
close all; clear all; clc;

%% Setup Parameters
% Step size in microns
step = 1;

% Origin
origin_index = 801;
%798
% B4C after origin_index
% SiC before origin_index

% Interval for moving average
interval = 51;

% Specify interval over which to search for peak maximum
sample_name = 'D4';
sample_material = 'B4C';
% sample_material = 'SiC';
    
if strcmp(sample_material, 'B4C')
    peak_id = {'480', '530', '730', '1080'};
    search_loc = [460,500; 510,550; 700,760; 1040,1120];
elseif strcmp(sample_material, 'SiC')
    peak_id = {'790', '970'};
    search_loc = [770,810; 950,990];
end

fig_text = cell(1,length(peak_id));
title_text = cell(1,length(peak_id));

for i = 1:length(peak_id)
    fig_text(1,i) = cellstr(append(sample_name, ' - ', sample_material, ' - ', peak_id{1,i}, ' peak.png'));
    title_text(1,i) = cellstr(append('\textbf{Boron Carbide - ', peak_id{1,i}, ' cm$^{-1}$ Peak Shift}'));
end

%% Data Selection Prompt
[file, path] = uigetfile({'*.txt'});
data_name = fullfile(path, file);
data = table2array(readtable(data_name));

if contains(path,'/')
    delim = '/';
elseif contains(path,'\')
    delim = '\';
end 

parent_folder = append(path, sample_name, delim);

mkdir([parent_folder, sample_material])

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

%%


%% Graph Labels
subtitle_text = '(Adjacent to Undoped Boron Carbide)';
x_text = '\textbf{Distance from Interface ($\mu$m)}';
y_text = '\textbf{Raman Shift (cm$^{-1}$)}';

for i = 1:length(peak_id)
    %% Determine index of start and end location
    title_name = title_text{1,i};
    fig_name = fig_text{1,i};
        
    [~,start_index] = min(abs(Wavenumber - search_loc(i,1)));
    [~,end_index] = min(abs(Wavenumber - search_loc(i,2)));

    disp_range = start_index:end_index;

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
    
    fig = figure('visible','off');
    subplot(2,1,1)
    plot(Distance, peak_loc(1,index), 'o')
    title(title_name, 'interpreter', 'latex', 'FontSize', 18)
    %subtitle(subtitle_text, 'interpreter', 'latex', 'FontSize', 14)
    xlabel(x_text, 'interpreter', 'latex', 'FontSize', 14)
    ylabel(y_text, 'interpreter', 'latex', 'FontSize', 14)
    xlim([Distance(1) Distance(end)])
    
    moving_avg = movmean(peak_loc(1,index),interval);
    
    subplot(2,1,2)
    plot(Distance, moving_avg, 'b', 'LineWidth', 2)
    title(title_name, 'interpreter', 'latex', 'FontSize', 18)
    %subtitle(subtitle_text, 'interpreter', 'latex', 'FontSize', 14)
    xlabel(x_text, 'interpreter', 'latex', 'FontSize', 14)
    ylabel(y_text, 'interpreter', 'latex', 'FontSize', 14)
    xlim([Distance(1) Distance(end)])

    new_path = append(path, delim, sample_name, delim, sample_material, delim);

    saveas(fig, fullfile(new_path, fig_name))
    
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
end