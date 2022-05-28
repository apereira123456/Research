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

%% Baseline Subtraction
l = length(Intensity);
lp = ceil(0.5*l);

initial_Spectrum = [ones(lp,1) * Intensity(1); Intensity; ones(lp,1) * Intensity(l)];

l2 = length(initial_Spectrum);
S = initial_Spectrum;
n = 1;
flag = 0;

while flag == 0
    n = n + 2;
    i = (n-1)/2;

    stripping = 0;

    y = sgolayfilt(S,0,n);

    m = length(S);

    Baseline = zeros(m,1);

    for i = 1:1:m
        if S(i) > y(i)
            stripping = 1;
            Baseline(i) = y(i);
        else
            Baseline(i) = S(i);
        end
    end

    A(i) = trapz(S-Baseline);

    Stripped_Spectrum{i} = Baseline;

    S = Baseline;

    if i > 3
        if A(i-1) < A(i-2) && A(i-1) < A(i)
            i_min = i-1;
            flag = 1;
        end
    end
end

Base = Stripped_Spectrum{i_min};
Corrected_Spectrum = initial_Spectrum - Base; Corrected_Spectrum = Corrected_Spectrum(lp+1:lp+l);
Base = Base(lp+1:lp+l);

Norm = normalize(Corrected_Spectrum,1,'range');

%     %% Baseline Subtraction 
%     [Cp,Sl,Ic] = ischange(Intensity,'linear');                          % Detect Changes, Calculates Slopes (& Intercepts)
%     [Cts,Edg,Bin] = histcounts(Sl, 10);                                 % Histogram Of Slopes
%     [Max,Binmax] = max(Cts);                                            % Find Largest Bin
%     LinearRegion = (Bin==Binmax);                                       % Logical Vector Of Values Corresponding To Largest Number Of Slopes
%     B = polyfit(Wavenumber(LinearRegion), Intensity(LinearRegion), 1);  % Linear Fit
%     L = polyval(B, Wavenumber);                                         % Evaluate
%     yc = Intensity - L; 
%     Norm = normalize(yc,1,'range');

%% 2D Plot
figure()
hold on
plot(Wavenumber, Norm, 'b', 'LineWidth', 2)
title(title_text, 'interpreter', 'latex', 'FontSize', 18)
xlabel(x_text, 'interpreter', 'latex', 'FontSize', 14)
ylabel(y_text, 'interpreter', 'latex', 'FontSize', 14)