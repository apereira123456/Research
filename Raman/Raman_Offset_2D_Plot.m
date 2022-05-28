%% Clear
close all; clear all; clc;

%% Setup Parameters
i = [0, 0];
j = [0, 648];
i_size = length(i);
peak = 1091;
range = [1050 1140];

%% Graph Labels
title_text = '\textbf{B$_4$C Spectrum}';
x_text = '\textbf{Raman Shift (cm$^{-1}$)}';
y_text = '\textbf{Intensity}';

%% User Data Selection Prompt
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
Wavenumber = zeros(s_size,1);
Intensity = zeros(s_size,1);

for u = 1:i_size
    for v = 1:s_size
        n = v + i(u)*s_size + j(u)*x_size*s_size;
        Wavenumber(v,u) = data(v,3);
        Intensity(v,u) = data(n,4);
        v = v + 1;
    end
end

%% Display Range
start_wavenumber = Wavenumber(s_size-1,1)
end_wavenumber = Wavenumber(1,1)

%% Display Range Prompt
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

Wavenumber = Wavenumber(disp_range,:);
Intensity = Intensity(disp_range,:);
Norm = normalize(Intensity,1,'range');

%% 2D Plot
figure()
subplot(2,1,1) %subplot(rows,coloumns,current_axis)
hold on
plot(Wavenumber(:,1), Norm(:,1), 'b', 'LineWidth', 2)
plot([peak,peak], [0,1], 'r', 'Linewidth', 2)

title('\textbf{B$_4$C 1080 cm$^{-1}$ Peak}', 'interpreter', 'latex', 'FontSize', 18)
subtitle('(Away from Interface)', 'interpreter', 'latex', 'FontSize', 14)
ylabel(y_text, 'interpreter', 'latex', 'FontSize', 14)
xlim(range)
ylim([0 1])

subplot(2,1,2)
hold on
plot(Wavenumber(:,2), Norm(:,2), 'b', 'LineWidth', 2)
plot([peak,peak], [0,1], 'r', 'Linewidth', 2)

title('\textbf{B$_4$C 1080 cm$^{-1}$ Peak}', 'interpreter', 'latex', 'FontSize', 18)
subtitle('(Near Interface)', 'interpreter', 'latex', 'FontSize', 14)
xlabel(x_text, 'interpreter', 'latex', 'FontSize', 14)
ylabel(y_text, 'interpreter', 'latex', 'FontSize', 14)
xlim(range)
ylim([0 1])