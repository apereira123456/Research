%% Clear
close all; clear all; clc;

%% Setup Parameters
% Specify excitation laser wavelength in nm
lambda_excitation = 532;

% Specify spectral bands to be analyzed
wave_start = 130;                   % spectral bands to be analyzed
wave_end = 1500;                    % spectral bands to be analyzed

%% Data Selection Prompt
[file, path] = uigetfile({'*.txt'});

UIFigure.Visible = 'on';
DataEditField.Value = file;

data_name = fullfile(path, file);
tic
data = table2array(readtable(data_name));
toc


%% Determine data array size
tic
[~,x,~] = unique(data(:,1),'rows');
[~,y,~] = unique(data(:,2),'rows');
[~,s,~] = unique(data(:,3),'rows');

[x_size, ~] = size(x);
[y_size, ~] = size(y);
[s_size, ~] = size(s);
toc

tic
s_start = data(1,3);
s_size = find(data(2:end,3) == s_start, 1, 'first');

x_size = (data(end,1) - data(1,1)) / (data(s_size + 1,1) - data(1,1)) + 1;
y_size = (data(end,2) - data(1,2)) / (data(x_size*s_size + 1,2) - data(1,2)) + 1;
toc

% %% Convert raw data to usable format
% % Store Raman shift
% wavenumber = data(s_size:-1:1,3);    
% 
% % Initialize the cell that stores the Raman data
% dcube = zeros(y_size,x_size,s_size);
% 
% tic
% for j = 1:y_size
%     for i = 1:x_size
%         a = (i-1)*s_size + (j-1)*x_size*s_size + 1;
%         b = i*s_size + (j-1)*x_size*s_size;
% 
%         dcube(j,i,:) = data(b:-1:a,4);
%     end
% end
% toc
% 
% tic
% avg = flipud(squeeze(mean(dcube, [1 2])));
% toc