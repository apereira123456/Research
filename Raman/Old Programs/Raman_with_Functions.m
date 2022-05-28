clearvars

%% User Data Selection Prompt
data_name = get_data();
data = table2array(readtable(data_name));

%% Logic for Scans of Different Dimensions
% Point Scan
if size(data,2) == 2
    Wavenumber = data(:,1);
    Intensity = data(:,2);
    
    figure()
    plot(Wavenumber, Intensity)

else
    [x_size, y_size, s_size] = scan_size(data);

    % Line Scan
    if x_size == 1 || y_size == 1
        i = zeros(s_size,1);
        j = (0:(s_size-1))';
    
    % Map Scan
    else
        image_name = get_image();
        image = imread(image_name);
        
        figure(), imshow(image)
        
        [x_coord, y_coord] = ginput(2);
        
        [i, j] = line_intersection(data, image, x_coord, y_coord);
    end
    
    i_size = length(i);
    X = zeros(s_size,1);
    Y = zeros(s_size,1);
    [Wavenumber, Intensity, Transition] = deal(zeros(s_size, i_size), zeros(s_size, i_size), zeros(s_size, i_size));
    
    % Convert txt file to plottable format
    for a = 1:i_size
            [X, Y] = text_to_curve(data, i(a), j(a));
            Wavenumber(:,a) = X;
            Intensity(:,a) = Y;
            Transition(:,a) = a;
    end

    
    %% 2D Plot of Scan
    figure()
    for a = 1:i_size
        hold on
        plot(Wavenumber(:,a), Intensity(:,a))
    end
    
    %% 3D Plot of Line Scan
    if x_size ~= 1 || y_size ~= 1
        % Determine scan range
        start_wavenumber = Wavenumber(s_size-1,1);
        end_wavenumber = Wavenumber(1,1);
        
        % Enter the range associated with the wavenumbers you would like displayed
        disp_range = 500:1000;
        
        %% Nonlinear Color Map
        % Set the bounds of the color map
        cMap = parula(256);
        [low, high] = bounds(data,1);
        min_intensity = low(4);
        max_intensity = high(4);
        most_often = mode(data,1);
        ctr_intensity = most_often(4);
        scaling_factor = 3;
        
        % Create color map with interpolated values
        x = 1:length(cMap);
        x = x - (ctr_intensity-min_intensity)*length(x)/(max_intensity-min_intensity);
        x = scaling_factor * x/max(abs(x));
        
        % Transform the color map so it is nonlinear 
        x = sign(x).* exp(abs(x));
        x = x - min(x); x = x*511/max(x)+1; 
        newMap = interp1(x, cMap, 1:512);
        
        %% 3D Mesh of Line Scan
        figure()
        colormap(newMap)
        surface = mesh(Wavenumber(disp_range,:), Transition(disp_range,:), Intensity(disp_range,:));
        surface.FaceColor = 'none';
        surface.EdgeColor = 'interp';
        xlabel('Wavelength (cm^-1)')
        ylabel('Transition')
        zlabel('Intensity')
        colorbar
    end
end