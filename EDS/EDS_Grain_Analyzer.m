%% Clear
close all; clear all; clc;

% Setup Parameters
go_again = 'Y';
n = 0;

% Atomic masses: B, C, Si (g/mol)
atomic_mass = [10.811, 12.011, 28.0855];

%% Specify scaling factors depending on SEM
SEM = 'Gemini';
% SEM = 'Sigma';

if strcmp(SEM,'Gemini')
    m = [0.75, -0.45];
    b = [0.83, 0.63];
elseif strcmp(SEM,'Sigma')
    m = [0.75, -0.45];
    b = [0.83, 0.63];
end

while go_again == 'Y'
    %% Specify sheetname
    prompt = "Enter sample name:";
    sheetname = input(prompt,'s');
    
    if isempty(sheetname)
        sheetname = 'Default';
    % check if sheetname exceeds Excel's 31 character limit and truncate it
    elseif size(sheetname,2)>31
        sheetname = sheetname(1:31);
    end
    
    %% Prompt user to select 3 data files
    % Start in default folder
    if n == 0
        start_path = '*.csv';
    end
    
    [filename, path] = uigetfile(start_path,'Select EDS data files','MultiSelect','on');
    number_of_files = length(filename);

    folders_in_path = regexp(path,'[\\/]','split');
    num_folders = length(folders_in_path);
    
    if contains(path,'/')
        delim = '/';
    elseif contains(path,'\')
        delim = '\';
    end 
    
    folders_in_path(num_folders-1:num_folders) = '';

    if n == 0     
        write_path = append(strjoin(folders_in_path,delim),delim);
        start_path = append(write_path,'*.csv');
    else
        start_path = append(strjoin(folders_in_path,delim),delim,'*.csv');
    end
    
    % Define filename to write data to
    excel_filename = 'EDS Grain Analysis.xlsx';
    full_excel_filename = fullfile(write_path, excel_filename);
    
    %% Read Data
    map_size = size(readmatrix(strcat(path,filename{1})));
    data = zeros(map_size(1),map_size(2),number_of_files);

    for i = 1:number_of_files
        % create path to file
        mapFile = strcat(path,filename{i});
    
        % read the data file
        data(:,:,i) = readmatrix(mapFile); % loads data from the file into a matrix
    end
    
    %% Grain selection prompt
    prompt = "Enter number of grains:";
    num_grains = input(prompt);
    
    [x_start, x_end, y_start, y_end] = deal(zeros(num_grains, 1));
    
    figure()
    hold on
    imshow(data)
    
    for i = 1:num_grains
        % Choose corners of rectangle
        [x_coord, y_coord] = ginput(2);
    
        x1 = x_coord(1); x2 = x_coord(2);
        y1 = y_coord(1); y2 = y_coord(2);
        
        % Determine corner and dimensions
        x = min(x1, x2);
        y = min(y1, y2);
        w = abs(x2 - x1);
        h = abs(y2 - y1);
        
        % Plot rectangle on top of image
        rectangle('Position', [x y w h], 'EdgeColor', 'r', 'LineWidth', 2)
    
        % Store values
        % Start is bottom left and end is top right
        x_start(i) = ceil(x); x_end(i) = floor(x + w);
        y_start(i) = ceil(y); y_end(i) = floor(y + h);
    end
    
    %% Compute average intensity for each grain
    [counts, intensity_ratio, scaling_factor, rel_intensity, wt_fraction, mol, at_percent] = deal(ones(num_grains,3));

    for i = 1:num_grains
        counts(i,1) = mean(data(y_start(i):y_end(i), x_start(i):x_end(i), 1), 'all');
        counts(i,2) = mean(data(y_start(i):y_end(i), x_start(i):x_end(i), 2), 'all');
        counts(i,3) = mean(data(y_start(i):y_end(i), x_start(i):x_end(i), 3), 'all');
    end

    %% Apply the Scaling Factors
    for i = 2:3
        % compute elemental intensity ratios to be used in scaling factor calcs
        intensity_ratio(:,i) = counts(:,i) ./ counts(:,1);
        % calculate scaling factors for each point along line averaged intensity
        scaling_factor(:,i) = b(i-1) * exp(m(i-1) .* intensity_ratio(:,i));
        % apply scaling factors to get corrected relative intensity of each element
        rel_intensity(:,i) = scaling_factor(:,i) .* intensity_ratio(:,i);
    end
    
    total_intensity = rel_intensity(:,1) + rel_intensity(:,2) + rel_intensity(:,3);
    
    for i = 1:3
        % convert corrected intensities to weight fractions
        wt_fraction(:,i) = rel_intensity(:,i) ./ total_intensity;
        % convert weight fractions to relative mol fraction
        mol(:,i) = wt_fraction(:,i) / atomic_mass(i);
    end
    
    % convert relative mol fraction to atomic percents
    total_mols = mol(:,1) + mol(:,2) + mol(:,3);
    
    for i = 1:3
        at_percent(:,i) = 100 * mol(:,i) ./ total_mols;
    end

    b_to_c_ratio = at_percent(:,1) ./ at_percent(:,2);
    
    %% Prepare data for export
    % combine the 2theta, counts, and normalized counts columns
    two_column_data = array2table(num2cell([at_percent(:,3), b_to_c_ratio]));
    
    % create data headers to label the columns]
    column_labels = {'Si (at. %)', 'B/C Ratio'};
    separator = {blanks(1), blanks(1)};
    
    % combine data and headers into one big table to write to Excel
    composite_data = [column_labels; two_column_data; separator];
    
    %% Write data to output Excel spreadsheet
    writetable(composite_data,full_excel_filename,'Sheet',sheetname,'WriteVariableNames',false,'WriteMode','append');
    
    % final printout
    fprintf('\nComplete. Sheets written to [%s] in [%s] \n',excel_filename,write_path)

    n = n + 1;

    prompt = "Would you like to analyze another scan (Y/N)";
    go_again = input(prompt,'s');

    if isempty(go_again)
        go_again = 'N';
    end
end