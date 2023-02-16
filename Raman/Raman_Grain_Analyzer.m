%% Clear
close all; clear all; clc;

%% Setup Parameters
% Specify excitation laser wavelength in nm
lambda_excitation = 532;

% Specify spectral bands to be analyzed
wave_start = 130;                   % spectral bands to be analyzed
wave_end = 1500;                    % spectral bands to be analyzed

% Specify RGB bands
band = [42, 121, 384];

num_g = 0;

%% Data Selection Prompt
[file, path] = uigetfile({'*.txt'});

UIFigure.Visible = 'on';
DataEditField.Value = file;

data_name = fullfile(path, file);
data = readmatrix(data_name);


%% Determine data array size
s_size = find(data(2:end,3) == data(1,3), 1, 'first');

x_size = round((data(end,1) - data(1,1)) / (data(s_size + 1,1) - data(1,1)) + 1);
y_size = round((data(end,2) - data(1,2)) / (data(x_size*s_size + 1,2) - data(1,2)) + 1);
ratio = x_size / y_size;

%% Convert raw data to usable format
% Store Raman shift
wavenumber = data(s_size:-1:1,3);    

% Initialize the cell that stores the Raman data
dcube = zeros(y_size,x_size,s_size);

for j = 1:y_size
    for i = 1:x_size
        a = (i-1)*s_size + (j-1)*x_size*s_size + 1;
        b = i*s_size + (j-1)*x_size*s_size;

        dcube(j,i,:) = data(b:-1:a,4);
    end
end

[~, start_index] = min(abs(wavenumber - wave_start));
[~, end_index] = min(abs(wavenumber - wave_end));

num_bands = end_index - start_index + 1;

wavenumber = wavenumber(start_index:end_index);
dcube = normalize(dcube(:,:,start_index:end_index),3,'range');

%% Grain Visualization
for i = 1:3
    grain_map(:,:,i) = histeq(dcube(:,:,band(i)));
end

figure
img = imshow(grain_map);

%% Polygon ROI
roi = images.roi.Polygon(gca);
draw(roi)

num_g = num_g + 1;

[pos] = roi.Position;
x = pos(:,1);
y = pos(:,2);
s = num2str(num_g);

pshape = polyshape(x,y);
[a,b] = centroid(pshape);
text(a,b,s, 'FontSize', 24, 'FontWeight', 'bold', 'Color', 'w')

grain_data{num_g} = mat2cell(pos, length(x), 2);

%%
roi = poly2mask(x,y,y_size,x_size);
    
temp_roi = false(y_size,x_size,num_bands);
temp_roi(:,:,band(2)) = roi;

Peak_480 = sum(dcube(temp_roi), 'all') / sum(roi, 'all');


%% Write Data
for i = 1:num_g
    x = grain_data{i}(:,1);
    y = grain_data{i}(:,2);
    roi = poly2mask(x,y,y_size,x_size);
    
    temp_roi = zeros(y_size,x_size,num_bands);
    temp_roi(:,:,band(2)) = roi;

    Grain(i) = i;
    Peak_480(i) = sum(dcube(temp_roi), 'all') / sum(roi, 'all');
    
    img_mask = repmat(roi, [1 1 3]);
    img = img_mask .* grain_map;
    
    s = num2str(i);
    figure_name = fullfile(path, ['GrainImages/Grain ' s '.png']);
    imwrite(img, figure_name);
end
    
%% Write Table
peak_intensity = table(Grain,Peak_480);

table_name = append(path, 'Grain_Data.txt');
writetable(peak_intensity, table_name, 'Delimiter', 'tab')