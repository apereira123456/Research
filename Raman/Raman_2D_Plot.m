%% Clear
close all; clear all; clc;

%% Setup Parameters
i = 0;
j = 798;

%% User Data Selection Prompt
[file, path] = uigetfile({'*.txt'});
data_name = fullfile(path, file);
data = table2array(readtable(data_name));

%% Logic for Scans of Different Dimensions
% Point Scan
if size(data,2) == 2
    start_wavenumber = data(end,1);
    end_wavenumber = data(1,1);

    Wavenumber = data(:,1);
    Intensity = data(:,2);

    fprintf('Start Wavenumber: %.2f \n',start_wavenumber)
    fprintf('End Wavenumber: %.2f \n',end_wavenumber)
    
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
    Wavenumber = zeros(s_size,1);
    Intensity = zeros(s_size,1);
    
    for v = 1:s_size
        n = v + i*s_size + j*x_size*s_size;
        Wavenumber(v,1) = data(v,3);
        Intensity(v,1) = data(n,4);
        v = v + 1;
    end
    
    %% Display Range Prompt
    start_wavenumber = Wavenumber(s_size-1,1);
    end_wavenumber = Wavenumber(1,1);

    fprintf('Start Wavenumber: %.2f \n',start_wavenumber)
    fprintf('End Wavenumber: %.2f \n',end_wavenumber)
end

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

Wavenumber = Wavenumber(disp_range,1);
Intensity = Intensity(disp_range,1);

%% Graph Labels
title_text = '\textbf{B$_4$C Spectrum}';
x_text = '\textbf{Raman Shift (cm$^{-1}$)}';
y_text = '\textbf{Intensity}';

Norm = normalize(Intensity,1,'range');

%% 2D Plot
figure()
hold on
plot(Wavenumber, Norm, 'b', 'LineWidth', 2)
title(title_text, 'interpreter', 'latex', 'FontSize', 18)
xlabel(x_text, 'interpreter', 'latex', 'FontSize', 14)
ylabel(y_text, 'interpreter', 'latex', 'FontSize', 14)