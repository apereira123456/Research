% Kent Christian 11/5/21
% Modified by Andrew Pereira 05/17/2022
% Read .xrdml files of XRD data and turn them into Excel spreadsheets
% Can select multiple .xrdml files to generate to 1 output file

close all; clear all; clc;

% ask user to locate data file then automatically create output file name
[filenames, path] = uigetfile('*','MultiSelect','on','Select one or more XRD data files (.xrdml)');
%xrd_file = strcat(path,filenames);
[~, filename_no_extension, ~] = fileparts(filenames);

% Determine name of folder the file is in
folders_in_path = regexp(path,'[\\/]','split');
size_path = size(folders_in_path,2);
folder_name = folders_in_path{size_path - 1};

% make output Excel named after datafile, or after folder if multiple
if isa(filenames,'cell')
    % multiple files selected. Filenames is cell
    filename = strcat(folder_name,'.xlsx');
    full_filename = strcat(fullfile(path, folder_name),'.xlsx');
else
    % only 1 file selected. Filenames is str
    filename = strcat(filename_no_extension,'.xlsx');
    full_filename = strcat(fullfile(path, filename_no_extension),'.xlsx');
    filenames = {filenames};
end

% iterate over each file the user selected to read and write to spreadsheet
for i = 1:size(filenames,2)
    % open and read the file
    raw_data = fileread(fullfile(path, filenames{i}));  

    % Extract the data from the file
    % find the tag that signifies where data starts in the file (start tag)
    start_tag = '<intensities unit="counts">';
    start_tag_length = size(start_tag,2);
    data_start_index = strfind(raw_data,start_tag) + start_tag_length;

    % find the tag that signifies where data ends in the file (end tag)
    end_tag = '</intensities>';
    data_end_index = strfind(raw_data,end_tag)-1;

    % grab just the intensity counts data, which is between start and end tags
    data_string = raw_data(data_start_index:data_end_index);
    counts = str2num(data_string)';
    normalized_counts = counts/max(counts);

    % Find scan range
    start_pos = '<positions axis="2Theta" unit="deg">';
    end_pos = '<positions axis="Omega" unit="deg">';
    angle_start_index = strfind(raw_data,start_pos) + 43;
    angle_end_index = strfind(raw_data,end_pos) - 25;

    range_string = raw_data(angle_start_index:angle_end_index);

    % Return only numbers from the scan range
    angle = str2double(regexp(range_string,'\d+\.?\d*','match'));

    % generate 2theta values automatically
    data_size = size(counts,1);
    two_theta = linspace(angle(1),angle(2),data_size)';

    % combine the 2theta, counts, and normalized counts columns
    three_column_data = array2table(num2cell([two_theta, counts, normalized_counts]));

    % create data headers to label the columns
    top_row = {filenames{i}(1:end-6), '', ''};
    column_labels = {'2theta (deg)', 'Intensity (counts)', 'Normalized Intensity (counts)'};
    header = array2table([top_row; column_labels]);

    % combine data and headers into one big table to write to Excel
    composite_data = [header; three_column_data];

    %% Write data to output Excel spreadsheet
    % check if sheetname exceeds Excel's 31 character limit and truncate it
    sheetname = filenames{i}(1:end-6);
    if size(sheetname,2)>31
        sheetname = filenames{i}(1:31);
    end

    writetable(composite_data,full_filename,'Sheet',sheetname,'WriteVariableNames',false);
    range_string = sprintf('Sheet %.f: %s',i,sheetname); disp(range_string)
end

% final printout
range_string = sprintf('\nComplete. Sheets written to [%s] in [%s]',filename,path); disp(range_string)
