%% Clear
close all; clear all; clc;

%% Setup Parameters
% Step size in microns
step = 1;

% Origin
origin_index = 800;

% Interval for moving average
interval = 51;

%
check = 1;

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

%% Convert raw data to plottable format
sorted_data = cell(x_size,y_size,2);
sorted_data(:,:,1) = mat2cell(data(s_size:-1:1,3),s_size,1); 

for j = 1:y_size
    for i = 1:x_size
        a = (i-1)*s_size + (j-1)*x_size*s_size + 1;
        b = i*s_size + (j-1)*x_size*s_size;

        sorted_data(i,j,2) = mat2cell(normalize(data(b:-1:a,4),1,'range'),s_size,1);
    end
end

% Line Scan 
% Can only parse vertical or horizontal lines
if x_size == 1 || y_size == 1
    %% Location Selection Prompt
    while check
        prompt = "Enter location start value:";
        loc_start = input(prompt);

        % Use the full line if no input is given
        if isempty(loc_start)
            loc_range = 1:(x_size + y_size - 2);
            check = 0;

        else
            prompt = "Enter location end value:";
            loc_end = input(prompt);
        
            loc_range = loc_start:loc_end;
            size(loc_range)

            % Check for valid range
            if loc_start >= 0 && loc_end <= x_size + y_size - 2 && loc_start <= loc_end
                check = 0;

            else
                disp('Invalid input.')
            end
        end
    end

    j = 1:y_size;
    i = ones(size(j,1));

%     i = 1:x_size;
%     j = ones(size(j,1));

%     i = zeros(length(loc_range),1);
%     j = loc_range';

% Map Scan
else
    
end

% X(:,u) = Wavenumber;
% Y(:,u) = Intensity;

%% Display range prompt
start_wavenumber = Wavenumber(1,1);
end_wavenumber = Wavenumber(end,1);

fprintf('Start Wavenumber: %.2f \n',start_wavenumber)
fprintf('End Wavenumber: %.2f \n',end_wavenumber)

check = 1;

while check
    prompt = "Enter wavenumber start value:";
    a = input(prompt);
    
    % Use the full range if no input is given
    if isempty(a)
        start_index = 1;
        end_index = s_size;
        disp_range = start_index:end_index;
        check = 0;
    
    % Otherwise find the index of the closest wavenumber to the input
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
[M, I] = max(Y(disp_range,:), [], 1);
peak_loc = X(disp_start + I);

if x_size == 1 || y_size == 1
    Distance = step * (loc_range - origin_index);
else
    Range = 0:length(loc_range) - 1;
    Distance = step * Range;    
end

figure()
plot(Distance, peak_loc, 'o')
title(title_text, 'interpreter', 'latex')
xlabel(x_text, 'interpreter', 'latex')
ylabel(y_text, 'interpreter', 'latex')

moving_avg = movmean(peak_loc,interval);
moving_std = movstd(peak_loc,interval);

top = moving_avg + moving_std;
bottom = moving_avg - moving_std;

figure()
hold on
plot(Distance, moving_avg, 'b', 'LineWidth', 2)
title(title_text, 'interpreter', 'latex', 'FontSize', 18)
%subtitle(subtitle_text, 'interpreter', 'latex', 'FontSize', 14)
xlabel(x_text, 'interpreter', 'latex', 'FontSize', 14)
ylabel(y_text, 'interpreter', 'latex', 'FontSize', 14)

x_std = [Distance, fliplr(Distance)];
y_std = [top, fliplr(bottom)];

% Standard Deviation
%     fill(x_std, y_std, [0.5 0.5 0.5], 'FaceAlpha', 0.3);