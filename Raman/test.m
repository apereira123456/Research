%%
load('fit_data.mat')

%%
% pixel headers
% M: material, X: x-coordinate, Y: y-coordinate, H: table height

% peak headers
% 'Position' 'Height' 'FWHM' 'Area'

%%
for i = 1:x_size
    y = 1;
    for j = 1:y_size
        index = i + (j-1) * x_size;

        if isempty(fit_data{j,i})
            height = 0;
            conversion(y:y + height, 1 + 4 * (i - 1):4 + 4 * (i - 1)) = [index height NaN NaN];
        else
            height = size(fit_data{j,i},1);
            conversion(y:y + height, 1 + 4 * (i - 1):4 + 4 * (i - 1)) = [index height NaN NaN; fit_data{j,i}(:,2:5)];
        end
        y = y + height + 1;
    end
end

conversion(1,end+1) = x_size;
conversion(end+1,1) = y_size;

writematrix(conversion,'conversion.txt','Delimiter','tab')

%%
table = table2array(readtable('conversion.txt'));

x_size = table(1,end);
y_size = table(end,1);

test_data = cell(y_size,x_size);

for i = 1:x_size
    y = 1;
    for j = 1:y_size
        height = table(y, 2 + 4 * (i - 1));

        test_data(j,i) = mat2cell(table(y+1:y + height, 1 + 4 * (i - 1):4 + 4 * (i - 1)), height, 4);
        
        y = y + height + 1;
    end
end

%%
peak_map = NaN(y_size,x_size);

for j = 1:y_size
    for i = 1:x_size
        if classMap(j,i) == 1
            peak_map(j,i) = test_data{j,i}(3,2);
        end
    end
end 

%%
varTypes = {'string', 'uint8', 'uint8', 'uint8', 'double', 'double', 'double', 'double'};
varNames = {'M', 'X', 'Y', 'H', 'Position', 'Height', 'FWHM', 'Area'};
peak_export = table('Size',[1 8],'VariableTypes',varTypes,'VariableNames',varNames);
peak_export{1,:} = 'B4C', 1, 1, 1, 1087.4, 1, 32.3, 16.15;