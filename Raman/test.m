%% Standard Selection Prompt
UIFigure.Visible = 'off';
[file, temp_path] = uigetfile({'*.txt'});

UIFigure.Visible = 'on';
PeakFitEditField.Value = file;

fit_data_name = fullfile(temp_path, file);            
peak_import = table2array(readtable(fit_data_name));

x_size = max(peak_import(:,2));
y_size = max(peak_import(:,3));

phase_map = zeros(y_size, x_size);
peak_fit = cell(y_size,x_size);

table_row = 1;

for j = 1:y_size
    for i = 1:x_size
        if peak_import(table_row,1) == 1
            phase_map(j,i) = peak_import(table_row,1);
            peak_fit(j,i) = mat2cell(peak_import(table_row:table_row + 6,5:8),7,4);
            table_row = table_row + 8;           
        
        elseif peak_import(table_row,1) == 4
            phase_map(j,i) = peak_import(table_row,1);
            peak_fit(j,i) = mat2cell(peak_import(table_row,5:8),1,4);
            table_row = table_row + 2;
            
        else
            phase_map(j,i) = peak_import(table_row,1);
            table_row = table_row + 2;
        end
    end
end

phase_map1 = phase_map;

% phase_map1(phase_map == 0) = 2 * (phase_map == 0);
% phase_map1(phase_map == 1) = 8 * (phase_map == 1);
% phase_map1(phase_map == 2) = 2 .* phase_map1(phase_map == 2);
% phase_map1(phase_map == 3) = 2 .* phase_map1(phase_map == 3);
% phase_map1(phase_map == 4) = 2 * (phase_map == 4);

b4c_num = sum(phase_map == 1,'all');
sic_num = sum(phase_map == 4,'all');
other_num = sum(phase_map == 0,'all');
total = x_size * y_size;

total(:,:,1) = sum(8 * (phase_map == 1),'all');
total(:,:,2) = sum(2 * (phase_map == 4),'all');
total(:,:,3) = sum(2 * (phase_map == 0),'all');
total_row = sum(total,'all');

[M, X, Y, H, Position, Height, FWHM, Area] = deal(zeros(sum(phase_map1,'all'),1));