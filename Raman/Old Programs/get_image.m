function [image_name] = get_image()
    [file, path] = uigetfile({'*.png;*.jpg;*.jpeg;*.tif'});
    image_name = fullfile(path, file);
end