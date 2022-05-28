function [data_name] = get_data()
    [file, path] = uigetfile({'*.txt'});
    data_name = fullfile(path, file);
end