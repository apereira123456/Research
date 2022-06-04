%% Clear
close all; clear all; clc;

%% Setup Parameters
% Specify index of data you wish to plot
index = 800;

%
check = 1;

%% User Data Selection Prompt
[file, path] = uigetfile({'*.txt'});
data_name = fullfile(path, file);
data = table2array(readtable(data_name));

%% Logic for Scans of Different Dimensions
% Point Scan
if size(data,2) == 2
    Wavenumber = data(:,1);
    Intensity = data(:,2);
  
% Line or Map Scan
else
    %% Find Scan Size
    [~,x,~] = unique(data(:,1),'rows');
    [~,y,~] = unique(data(:,2),'rows');
    [~,s,~] = unique(data(:,3),'rows');
    
    [x_size, ~] = size(x);
    [y_size, ~] = size(y);
    [s_size, ~] = size(s);
    
    %% Extract Plot from Data
    i = mod(index, x_size);
    if i == 0
        i = x_size;
    end
    j = ceil(index/x_size);

    a = (i-1)*s_size + (j-1)*x_size*s_size + 1;
    b = i*s_size + (j-1)*x_size*s_size;

    Wavenumber = data(s_size:-1:1,3);
    Intensity = normalize(data(b:-1:a,4),1,'range');
end

%% Display Range Prompt
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
title_text = '\textbf{B$_4$C Spectrum}';
x_text = '\textbf{Raman Shift (cm$^{-1}$)}';
y_text = '\textbf{Intensity}';

%% 2D Plot
figure()
plot(Wavenumber(disp_range), Intensity(disp_range,1), 'b', 'LineWidth', 2)

title(title_text, 'interpreter', 'latex', 'FontSize', 18)
xlabel(x_text, 'interpreter', 'latex', 'FontSize', 14)
ylabel(y_text, 'interpreter', 'latex', 'FontSize', 14)
xlim([Wavenumber(start_index,1) Wavenumber(end_index,1)])
ylim([0 1])