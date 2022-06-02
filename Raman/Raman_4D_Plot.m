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

%% Determine data array size
[~,x,~] = unique(data(:,1),'rows');
[~,y,~] = unique(data(:,2),'rows');
[~,s,~] = unique(data(:,3),'rows');

[x_size, ~] = size(x);
[y_size, ~] = size(y);
[s_size, ~] = size(s);

%% Convert raw data to plottable format
i_size = length(i);

Wavenumber = zeros(s_size,1);
Intensity = zeros(s_size,1);
[X, Y, Z] = deal(zeros(s_size, i_size), zeros(s_size, i_size), zeros(s_size, i_size));

if x_size == 1 || y_size == 1
    Distance = step * (loc_range - origin_index);
else
    Range = 0:i_size - 1;
    Distance = step * Range;    
end

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

for u = 1:s_size
    Z(u,:) = Distance;
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
title_text = '\textbf{SiC - 1.5 at.\% Si + B$_{6.5}$C Diffusion Couple}';
x_text = '\textbf{Raman Shift (cm$^{-1}$)}';
y_text = '\textbf{Distance from Interface ($\mu$m)}';
z_text = '\textbf{Intensity}';

%% Normalized 3D Mesh of Line Scan
Norm = normalize(Y(disp_range,:),1,'range');

% Set the bounds of the color map
cMap = parula(256);
[low, high] = bounds(Norm,1);
min_intensity = low;
max_intensity = high;
most_often = mode(Norm,1);
ctr_intensity = most_often;
scaling_factor = 2;

% Create color map with interpolated values
x = 1:length(cMap);
x = x - (ctr_intensity - min_intensity)*length(x)/(max_intensity - min_intensity);
x = scaling_factor * x/max(abs(x));

% Transform the color map so it is nonlinear 
x = sign(x).* exp(abs(x));
x = x - min(x); 
x = x*511/max(x)+1; 
norm_map = interp1(x, cMap, 1:512);

% Colored faces with white edges
figure()
colormap(norm_map)
surface = mesh(X(disp_range,:), Z(disp_range,:), Norm);
surface.FaceColor = 'interp';
surface.EdgeColor = 'white';
surface.EdgeAlpha = '0.3';
title(title_text, 'interpreter', 'latex')
xlabel(x_text, 'interpreter', 'latex')
ylabel(y_text, 'interpreter', 'latex')
zlabel(z_text, 'interpreter', 'latex')