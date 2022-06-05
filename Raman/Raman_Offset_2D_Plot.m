%% Clear
close all; clear all; clc;

%% Setup Parameters
% Specify indices of data you wish to plot
index = [802, 1601];

% Specify peak position
peak = 1091;

%
check = 1;

%% Data Selection Prompt
[file, path] = uigetfile({'*.txt'});
data_name = fullfile(path, file);
data = table2array(readtable(data_name));

%% Find Scan Size
[~,x,~] = unique(data(:,1),'rows');
[~,y,~] = unique(data(:,2),'rows');
[~,s,~] = unique(data(:,3),'rows');

[x_size, ~] = size(x);
[y_size, ~] = size(y);
[s_size, ~] = size(s);

%% Extract Plot from Data
Wavenumber = zeros(s_size,2);
Intensity = zeros(s_size,2);

%% Extract Plot from Data
for n = 1:2
    i = mod(index(n), x_size);
    if i == 0
        i = x_size;
    end
    j = ceil(index(n)/x_size);
    
    a = (i-1)*s_size + (j-1)*x_size*s_size + 1;
    b = i*s_size + (j-1)*x_size*s_size;
    
    Wavenumber(:,n) = data(s_size:-1:1,3);
    Intensity(:,n) = normalize(data(b:-1:a,4),1,'range');
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
        [~,start_index] = min(abs(Wavenumber(:,1) - a));
        prompt = "Enter wavenumber end value:";
        a = input(prompt);
        [~,end_index] = min(abs(Wavenumber(:,1) - a));
    
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
subplot(2,1,1) %subplot(rows,coloumns,current_axis)
hold on
plot(Wavenumber(disp_range,1), Intensity(disp_range,1), 'b', 'LineWidth', 2)
plot([peak,peak], [0,1], 'r', 'Linewidth', 2)

title('\textbf{B$_4$C 1080 cm$^{-1}$ Peak}', 'interpreter', 'latex', 'FontSize', 18)
subtitle('(Away from Interface)', 'interpreter', 'latex', 'FontSize', 14)
ylabel(y_text, 'interpreter', 'latex', 'FontSize', 14)
xlim([Wavenumber(start_index,1) Wavenumber(end_index,1)])
ylim([0 1])

subplot(2,1,2)
hold on
plot(Wavenumber(disp_range,2), Intensity(disp_range,2), 'b', 'LineWidth', 2)
plot([peak,peak], [0,1], 'r', 'Linewidth', 2)

title('\textbf{B$_4$C 1080 cm$^{-1}$ Peak}', 'interpreter', 'latex', 'FontSize', 18)
subtitle('(Near Interface)', 'interpreter', 'latex', 'FontSize', 14)
xlabel(x_text, 'interpreter', 'latex', 'FontSize', 14)
ylabel(y_text, 'interpreter', 'latex', 'FontSize', 14)
xlim([Wavenumber(start_index,1) Wavenumber(end_index,1)])
ylim([0 1])