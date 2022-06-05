%% Clear
close all; clear all; clc;

%% Setup Parameters
% Step size in microns
step = 5;

% Origin
origin_index = 800;

%
check = 1;

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
[Wavenumber, Intensity] = deal(zeros(s_size,ij));

for n = index
    i = mod(n, x_size);
    if i == 0
        i = x_size;
    end
    j = ceil(n/x_size);

    a = (i-1)*s_size + (j-1)*x_size*s_size + 1;
    b = i*s_size + (j-1)*x_size*s_size;
 
    Wavenumber(:,n) = data(s_size:-1:1,3);
    Intensity(:,n) = normalize(data(b:-1:a,4),1,'range');
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

Distance = repelem(Distance,s_size,1);

%% Display Range Prompt
start_wavenumber = Wavenumber(1,index(1));
end_wavenumber = Wavenumber(end,index(1));

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
        [~,start_index] = min(abs(Wavenumber(:,index(1)) - a));
        prompt = "Enter wavenumber end value:";
        a = input(prompt);
        [~,end_index] = min(abs(Wavenumber(:,index(1)) - a));
    
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
title_text = '\textbf{SiC - 1.5 at.\% Si + B$_{6.5}$C Diffusion Couple}';
x_text = '\textbf{Raman Shift (cm$^{-1}$)}';
y_text = '\textbf{Distance from Interface ($\mu$m)}';
z_text = '\textbf{Intensity}';

%% Nonlinear Color Map Setup
% Set the bounds of the color map
cMap = parula(256);
[low, high] = bounds(Intensity(disp_range,index(1)),1);
min_intensity = low;
max_intensity = high;
most_often = mode(Intensity(disp_range,index(1)),1);
ctr_intensity = most_often;
scaling_factor = 2;

% Create color map with interpolated values
c = 1:length(cMap);
% Below matrix is rank deficient
c = c - (ctr_intensity - min_intensity)*length(c)/(max_intensity - min_intensity);
c = scaling_factor * c/max(abs(c));

% Transform the color map so it is nonlinear 
c = sign(c).* exp(abs(c));
c = c - min(c); 
c = c*511/max(c)+1; 
norm_map = interp1(c, cMap, 1:512);

%% 3D Mesh of Line Scan
% Colored faces with white edges
figure()
colormap(norm_map)
surface = mesh(Wavenumber(disp_range,index), Distance(disp_range,:), Intensity(disp_range,index));
surface.FaceColor = 'interp';
surface.EdgeColor = 'white';
surface.EdgeAlpha = '0.01';
title(title_text, 'interpreter', 'latex')
xlabel(x_text, 'interpreter', 'latex')
ylabel(y_text, 'interpreter', 'latex')
zlabel(z_text, 'interpreter', 'latex')