%% Clear
close all; clear all; clc;

%% Setup Parameters
% Step size in microns
step = 1;

% Origin
origin_index = 800;

%% User Data Selection Prompt
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

% Line Scan 
% Can only parse vertical or horizontal lines
if x_size == 1 || y_size == 1
    %% Location Selection Prompt
    check = 1;
    
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

    i = zeros(length(loc_range),1);
    j = loc_range';

% Map Scan
else
    %% Image selection prompt
    [file, path] = uigetfile({'*.png;*.jpg;*.jpeg;*.tif'});
    image_name = fullfile(path, file);
    image = imread(image_name);
    
    %% Coordinate selection prompt
    figure(), imshow(image)
    [x_coord, y_coord] = ginput(2);

    % Line Intersections
    % Determine image size and generate scale factor
    [y_pixels, x_pixels, ~] = size(image);
    scale = [x_pixels/(x_size-1), y_pixels/(y_size-1)];
    
    % Scale x-values and determine min/max
    x_pos = x_coord / scale(1);
    x_min = min(x_pos, [], 1);
    x_max = max(x_pos, [], 1);
    
    % Scale y-values and determine min/max
    y_coord = y_pixels - y_coord;
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
    ex_bdry = find(hv_int(:,1) < floor(x_min) | hv_int(:,1) < 0);
    ex_bdry = [find(hv_int(:,1) > ceil(x_max) | hv_int(:,1) > x_size - 1); ex_bdry];
    ex_bdry = [find(hv_int(:,2) < y_min | hv_int(:,2) < 0 ); ex_bdry];
    ex_bdry = [find(hv_int(:,2) > y_max | hv_int(:,2) > y_size - 1); ex_bdry];
    hv_int(ex_bdry,:) = [];
    
    % Remove duplicate values
    if m < 0
        unq = unique(hv_int,'rows');
        srtd = flip(sortrows(unq,[1 2],{'descend' 'ascend'}));
    else
        srtd = unique(hv_int,'rows');
    end

    i = srtd(:,1);
    j = srtd(:,2);
    a = [x_pos(1), x_pos(2)];
    b = [y_pos(1), y_pos(2)];

    % Draw line and grid
    figure()
    plot(a,b)
    xlim([0 x_size-1])
    ylim([0 y_size-1])
    set(gca,'xtick',0:1:x_size-1)
    set(gca,'ytick',0:1:y_size-1)
    grid on
end

%% Convert raw data to plottable format
i_size = length(i);

Wavenumber = zeros(s_size,1);
Intensity = zeros(s_size,1);
[X, Y, Z] = deal(zeros(s_size, i_size), zeros(s_size, i_size), zeros(s_size, i_size));

for u = 1:i_size
    v = 1;
    while v <= s_size
        n = v + i(u)*s_size + j(u)*x_size*s_size;
        Wavenumber(v,1) = data(v,3);
        Intensity(v,1) = data(n,4);
        v = v + 1;
    end
    X(:,u) = Wavenumber;
    Y(:,u) = Intensity;
end

%% Display range prompt
start_wavenumber = Wavenumber(s_size-1,1);
end_wavenumber = Wavenumber(1,1);

fprintf('Start Wavenumber: %.2f \n',start_wavenumber)
fprintf('End Wavenumber: %.2f \n',end_wavenumber)

check = 1;

while check
    prompt = "Enter wavenumber start value:";
    a = input(prompt);
    
    % Use the full range if no input is given
    if isempty(a)
        disp_range = 1:s_size;
        disp_start = 0;
        check = 0;

    else
        [~,disp_end] = min(abs(Wavenumber - a));
        prompt = "Enter wavenumber end value:";
        a = input(prompt);
        [~,disp_start] = min(abs(Wavenumber - a));
    
        disp_range = disp_start:disp_end;
        
        % Check for valid input
        if disp_start > disp_end
            disp('Invalid input.')

        elseif length(disp_range) > 1
            check = 0;

        else
            disp('Invalid input.')
        end
    end
end

%% Graph Labels
title_text = '\textbf{Silicon Carbide - 960 cm$^{-1}$ Peak Shift}';
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

moving_avg = movmean(peak_loc,151);
moving_std = movstd(peak_loc,151);

top = moving_avg + moving_std;
bottom = moving_avg - moving_std;

x_std = [Distance, fliplr(Distance)];
y_std = [top, fliplr(bottom)];

figure()
hold on
plot(Distance, moving_avg, 'b', 'LineWidth', 2)
title(title_text, 'interpreter', 'latex', 'FontSize', 18)
%subtitle(subtitle_text, 'interpreter', 'latex', 'FontSize', 14)
xlabel(x_text, 'interpreter', 'latex', 'FontSize', 14)
ylabel(y_text, 'interpreter', 'latex', 'FontSize', 14)

% Standard Deviation
%     fill(x_std, y_std, [0.5 0.5 0.5], 'FaceAlpha', 0.3);