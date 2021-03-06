%% Clear
close all; clear all; clc;

%% Setup Parameters
% Step size in microns
step = 2;

% Origin
origin_index = 402;
%401 404
% B4C after origin_index
% SiC before origin_index

% Interval for moving average
interval = 25;

trunc = [150, 1500];

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

if strcmp(sample_material, 'B4C')
    for i = 1:length(peak_id)
        fig_text(1,i) = cellstr(append(sample_name, ' - ', sample_material, ' - ', peak_id{1,i}, ' peak.png'));
        title_text(1,i) = cellstr(append('\textbf{Boron Carbide - ', peak_id{1,i}, ' cm$^{-1}$ Peak Shift}'));
    end
elseif strcmp(sample_material, 'SiC')
    for i = 1:length(peak_id)
        fig_text(1,i) = cellstr(append(sample_name, ' - ', sample_material, ' - ', peak_id{1,i}, ' peak.png'));
        title_text(1,i) = cellstr(append('\textbf{Silicon Carbide - ', peak_id{1,i}, ' cm$^{-1}$ Peak Shift}'));
    end
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

%% Convert raw data to plottable format
raw_data = cell(x_size,y_size,2);
raw_data(:,:,1) = mat2cell(data(s_size:-1:1,3),s_size,1); 

for j = 1:y_size
    for i = 1:x_size
        a = (i-1)*s_size + (j-1)*x_size*s_size + 1;
        b = i*s_size + (j-1)*x_size*s_size;

        raw_data(i,j,2) = mat2cell(data(b:-1:a,4),s_size,1);
    end
end

[~,trunc_start] = min(abs(raw_data{1,1,1} - trunc(1)));
[~,trunc_end] = min(abs(raw_data{1,1,1} - trunc(2)));
trunc_range = trunc_start:trunc_end;

trunc_data = cell(x_size,y_size,2); 

for j = 1:y_size
    for i = 1:x_size
        temp_wavenumber = raw_data{i,j,1};
        temp_intensity = raw_data{i,j,2};

        temp_wavenumber = temp_wavenumber(trunc_start:trunc_end,1);
        temp_intensity = temp_intensity(trunc_start:trunc_end,1);
        
        trunc_data(i,j,1) = mat2cell(temp_wavenumber,length(trunc_range),1);
        trunc_data(i,j,2) = mat2cell(normalize(temp_intensity,1,'range'),length(trunc_range),1);
    end
end

%% Specify analysis region
% Prompt user to select csv file from Raman Map to Line program
prompt = "Enter start position:";
a = input(prompt);

prompt = "Enter end position:";
b = input(prompt);

if a < b
    index = a:b;
else
    index = b:a;
end

Wavenumber = trunc_data{1,1,1};

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
    Peak_Loc = zeros(1,length(index));
    random_var = zeros(1,11);
    
    for v = index 
        for u = 1:11
            [M, I] = max(trunc_data{u,v,2}(disp_range,1), [], 1);
    
            random_var(1,u) = trunc_data{u,v,1}(start_index + I - 1,1);
        end
        Peak_Loc(1,v) = mean(random_var,2);
    end
    
    % Line Scans: Adjust distance based on origin
    Distance = step * (index - origin_index);
    
    fig = figure('visible','off');
    subplot(2,1,1)
    plot(Distance, Peak_Loc(1,index), 'o')
    title(title_name, 'interpreter', 'latex', 'FontSize', 18)
    %subtitle(subtitle_text, 'interpreter', 'latex', 'FontSize', 14)
    xlabel(x_text, 'interpreter', 'latex', 'FontSize', 14)
    ylabel(y_text, 'interpreter', 'latex', 'FontSize', 14)
    xlim([Distance(1) Distance(end)])
    
    moving_avg = movmean(Peak_Loc(1,index),interval);
    
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