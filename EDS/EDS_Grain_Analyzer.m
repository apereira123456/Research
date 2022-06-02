close all; clear all; clc;

% Setup Parameters
go_again = 'Y';
n = 0;
same_sample = 0;

%% Specify scaling factors depending on SEM
SEM = 'Gemini';

if strcmp(SEM,'Gemini')
    cb_a = 0.83;    cb_m = 0.75;
    sb_a = 0.63;    sb_m = -0.45;
elseif strcmp(SEM,'Sigma')
    cb_a = 0.83;    cb_m = 0.75;
    sb_a = 0.63;    sb_m = -0.45;
else
    fprintf('You must specify the SEM')
    go_again = 'N';
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

    % Check if the sample is the same
    if n == 0
        last_sheetname = sheetname;
    elseif strcmp(sheetname,last_sheetname)
        same_sample = 1;
    else
        last_sheetname = sheetname;
        same_sample = 0;
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
    elseif contains(path,'\\')
        delim = '\\';
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
    
    [data_y_size, data_x_size, ~] = size(data);
    
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
        x_start(i) = floor(x); x_end(i) = floor(x + w);
        y_start(i) = floor(y); y_end(i) = floor(y + h);
    end
    
    %% Compute average intensity for each grain
    if same_sample
        sheet_start = length(at_percent_si);
    else
        sheet_start = 0;
        [region_b_counts, region_c_counts, region_si_counts] = deal(zeros(num_grains,1));
    end

    for i = 1:num_grains
        region_b_counts(sheet_start + i) = mean(data(y_start(i):y_end(i), x_start(i):x_end(i), 1), 'all');
        region_c_counts(sheet_start + i) = mean(data(y_start(i):y_end(i), x_start(i):x_end(i), 2), 'all');
        region_si_counts(sheet_start + i) = mean(data(y_start(i):y_end(i), x_start(i):x_end(i), 3), 'all');
    end
    
    %% Apply the Scaling Factors
    
    % compute elemental intensity ratios to be used in scaling factor calcs
    cb_intensity_ratio = region_c_counts ./ region_b_counts;
    sib_intensity_ratio = region_si_counts ./ region_b_counts;
    
    % calculate scaling factors for each point along line averaged intensity
    cb_scaling_factor = 0.83*exp(0.75 .* cb_intensity_ratio);
    sb_scaling_factor = 0.63*exp(-0.45 .* sib_intensity_ratio);
    
    % apply scaling factors to get corrected relative intensity of each element
    rel_intensity_c = cb_scaling_factor .* cb_intensity_ratio;
    rel_intensity_si = sb_scaling_factor .* sib_intensity_ratio;
    rel_intensity_b = ones(1,size(region_b_counts,2)); % B is 1 as the reference
    
    % convert corrected intensities to weight fractions
    total_intensity = rel_intensity_c + rel_intensity_si + rel_intensity_b;
    wt_fraction_c = rel_intensity_c ./ total_intensity;
    wt_fraction_si = rel_intensity_si ./ total_intensity;
    wt_fraction_b = rel_intensity_b ./ total_intensity;
    
    % atomic masses
    c_mass = 12.011; % grams/mol
    si_mass = 28.0855; % grams/mol
    b_mass = 10.811; % grams/mol
    
    % convert weight fractions to relative mol fraction
    mol_c = wt_fraction_c / c_mass;
    mol_si = wt_fraction_si / si_mass;
    mol_b = wt_fraction_b / b_mass;
    
    % convert relative mol fraction to atomic percents
    total_mols = mol_c + mol_si + mol_b;
    at_percent_c = 100 * mol_c ./ total_mols;
    at_percent_si = 100 * mol_si ./ total_mols;
    at_percent_b = 100 * mol_b ./ total_mols;
    b_to_c_ratio = at_percent_b ./ at_percent_c;
    
    %% Prepare data for export
    % combine the 2theta, counts, and normalized counts columns
    two_column_data = array2table(num2cell([at_percent_si, b_to_c_ratio]));
    
    % create data headers to label the columns]
    column_labels = {'Si (at. %)', 'B/C Ratio'};
    
    % combine data and headers into one big table to write to Excel
    composite_data = [column_labels; two_column_data];
    
    %% Write data to output Excel spreadsheet
    writetable(composite_data,full_excel_filename,'Sheet',sheetname,'WriteVariableNames',false);
    
    % final printout
    fprintf('\nComplete. Sheets written to [%s] in [%s] \n',excel_filename,write_path)

    n = n + 1;

    prompt = "Would you like to analyze another scan (Y/N)";
    go_again = input(prompt,'s');

    if isempty(go_again)
        go_again = 'N';
    end
end