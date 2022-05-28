% Kent Christian 5/17/21
% Turn Raman spectrum data of a line map into a spreadsheet array file


% Guide: Smooth Raman line scan in WiRE. Save as .txt then run this

close all; clear all; clc;

% ask user to locate .txt data file
% then automatically create output file name
[filename, path] = uigetfile('*','Select Raman spectrum map file (.txt)');
spectraFile = strcat(path,filename);
[~, filename_no_extension, ~] = fileparts(filename);
output_excel_filename = strcat(filename_no_extension,'.xlsx');

% Initial data 
total_data = dlmread(spectraFile,'',1,0); %loads all data into x by 4 matrix
xrange = total_data(:,1); %list of all x positions
yrange = total_data(:,2); %list of all y positions

% Normally, data has 4 columns. If there is height data (Z) there are 5.
% Either (x,y,w,i) or (x,y,z,w,i)
totwavenumrange = total_data(:,end-1); %wavenumrange, repeated m by n times
intensitydata = total_data(:,end); %intensity values of Raman spectra

% find map size. Note either x or y should be 1 for vertical/horizontal line scan
x = size(unique(total_data(:,1)),1); % total number of step positions (pixels) in x direction
y = size(unique(total_data(:,2)),1); % total number of step positions (pixels) in y direction
data.boxsize=[x y]; % pixel dimensions of your image [m n]
numPoints = max(x,y); % how many points are in the line scan

% get list of just coordinates
coords = [xrange,yrange];

% Find wavenumber range by reading from start until first coordinate change
for n=1:size(coords,1)
    if not(isequal(coords(n,:), coords(n+1,:)))
        data.numwavenum=n;
        break
    end
end

% reorganize single column intensity data 
% move each coordinate's intensities into their own column
% 3 arrays will be output. Regular, normalized, and offset.
outputArray = zeros(n,numPoints+1); % initialize array same size as output data
outputArray(:,1) = totwavenumrange(1:n); % put wavenumbers in very first column
normalizedArray = outputArray; % a normalized array 
offsetArray = outputArray; % a normalized array with offsets for plotting
offsetAmount = 0.2;

for p = 1:numPoints
    % generate output array, normalized array, and offset array
    startpoint = n*p-n+1;
    endpoint = p*n;
    columnData = intensitydata(startpoint:endpoint);
    normColumnData = columnData/max(columnData);
    outputArray(:,p+1) = columnData;
    normalizedArray(:,p+1) = normColumnData;
    offsetArray(:,p+1) = normColumnData +(p-1)*offsetAmount;
end
fprintf('Calculation complete.\n');

%% Write regular data to Sheet 1
outputTable = array2table(outputArray);
writetable(outputTable,[path output_excel_filename],'Sheet','Raw Data','WriteVariableNames',false);

% Write normalized data on Sheet 2
normalizedTable = array2table(normalizedArray);
writetable(normalizedTable,[path output_excel_filename],'Sheet','Normalized','WriteVariableNames',false);

% Write offset data on Sheet 3 for plotting
offsetTable = array2table(offsetArray);
writetable(offsetTable,[path output_excel_filename],'Sheet','Offset values for plotting','WriteVariableNames',false);

% printout
fprintf('Written to Excel sheet: %s.\n',filename);
fprintf('In Directory: %s.\n',pwd);
fprintf('Sheet 1: Raw data values.\n');
fprintf('Sheet 2: Normalized values.\n');
fprintf('Sheet 3: Offset values.\n');
fprintf('\n');
