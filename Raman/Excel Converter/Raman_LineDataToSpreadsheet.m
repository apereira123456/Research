% Kent Christian 5/17/21
% Turn Raman spectrum data of a line map into a spreadsheet array file


% Guide: Smooth and normalize Raman map in WiRE. Save as .txt then run this
% program on it.

close all; clear all; clc;

% ask user to locate data file
% then automatically create output file name
[filename, path] = uigetfile('*','Select Raman spectrum map file (.txt)');
spectraFile = strcat(path,filename);
[~, filename_no_extension, ~] = fileparts(filename);
output_excel_filename = strcat(filename_no_extension,'.xlsx');

% Initial data 
total_data = dlmread(spectraFile,'',1,0); %loads all data into x by 4 matrix
xrange = total_data(:,1); %x positions
yrange = total_data(:,2); %y positions
totwavenumrange = total_data(:,3); %wavenumrange, repeated m by n times
intensitydata = total_data(:,4); %intensity values of Raman spectra
% Note: If Z data exists, it will be column 3. And wavenumber 4, intensity 5 

% find map size. Note either x or y should be 1 for vertical/horizontal line scan
x = size(unique(total_data(:,1)),1); % total number of step positions (pixels) in x direction
y = size(unique(total_data(:,2)),1); % total number of step positions (pixels) in y direction
data.boxsize=[x y]; % pixel dimensions of your image [m n]
numPoints = max(x,y); % how many points are in the line scan

% Divide your data based on the size of the wavenumber range
% into m by n, where m = # scan and n = # wavenumber
% Find the number of wavenumbers (x-axis) values in the selected range
for n=1:size(xrange,1)
    if xrange(n)<xrange(n+1)
        data.numwavenum=n;
        break
    end
end

% split intensity data (column 4) into a column for each point
outputArray = zeros(n,numPoints+1);
outputArray(:,1) = totwavenumrange(1:n); % write wavenumbers to first column
offsetArray = outputArray;
offsetAmount = 0.2;
for p = 1:numPoints
    startpoint = n*p-n+1;
    endpoint = p*n;
    outputArray(:,p+1) = intensitydata(startpoint:endpoint);
    offsetArray(:,p+1) = intensitydata(startpoint:endpoint)+(p-1)*offsetAmount;
end

%% Write intensity map data to  excel spreadsheet
outputTable = array2table(outputArray);
writetable(outputTable,[path output_excel_filename],'Sheet',filename_no_extension,'WriteVariableNames',false);
fprintf('\n%s\nwritten to:\n%s\n',filename,pwd);

% Write offset data on Sheet 2 for plotting
offsetTable = array2table(offsetArray);
writetable(offsetTable,[path output_excel_filename],'Sheet','Offset values for plotting','WriteVariableNames',false);
fprintf('Offset values written in Sheet 2\n');