function Raman_Interactive_2D_Plot
    %% Clear
    close all; clear all; clc;

    %% User Data Selection Prompt
    [file, path] = uigetfile({'*.txt'});
    data_name = fullfile(path, file);
    data = table2array(readtable(data_name));
    
    %% Logic for Scans of Different Dimensions
    % Determine data array size
    [~,x,~] = unique(data(:,1),'rows');
    [~,y,~] = unique(data(:,2),'rows');
    [~,s,~] = unique(data(:,3),'rows');
    
    [x_size, ~] = size(x);
    [y_size, ~] = size(y);
    [s_size, ~] = size(s);
    
    %% Convert raw data to plottable format
    sorted_data = cell(x_size,y_size,2);
    sorted_data(:,:,1) = mat2cell(data(1:s_size,3),s_size,1); 
    
    for j = 1:y_size
        for i = 1:x_size
            a = (i-1)*s_size + (j-1)*x_size*s_size + 1;
            b = i*s_size + (j-1)*x_size*s_size;
    
            sorted_data(i,j,2) = mat2cell(data(a:b,4),s_size,1);
        end
    end
    
    %% Graph Labels
    title_text = '\textbf{B$_4$C Spectrum}';
    x_text = '\textbf{Raman Shift (cm$^{-1}$)}';
    y_text = '\textbf{Intensity}';
    
    %% Extract I and J from index
    ij = x_size * y_size;
    index = 1;
    i = mod(index, x_size);
    if i == 0
        i = x_size;
    end
    j = ceil(index/x_size);
    
    %% Start Plot
    Wavenumber = sorted_data{1,1,1};
    Intensity = sorted_data{1,1,2};
    Norm = normalize(Intensity,1,'range');
    
    %% 2D Plot
    set(0,'units','pixels');
    screen_size = get(0,'ScreenSize');
    
    x = (screen_size(3) - screen_size(1))/6;
    y = (screen_size(4) - screen_size(2))/6;
    w = 4*x;
    h = 4*y;
    
    fig = uifigure;
    fig.Position = [x y w h];
    fig.Resize = 0;

    txt = uitextarea(fig);
    txt.Position = [0.9*w 0.5*h 0.06*w 0.05*h];
    txt.Value = {'Index=1', 'I=1; J=1'};

    sld = uislider(fig, 'ValueChangingFcn', @(sld,event) sliderMoving(event,sld));
    sld.Position = [0.1*w 0.1*h 0.8*w 0.1*h];
    sld.Value = 1;
    sld.Limits = [1 ij];
    
    btn_up = uibutton(fig, 'ButtonPushedFcn', @one_up);
    btn_up.Position = [0.02*w 0.55*h 0.07*w 0.05*h];
    btn_up.Text = 'One Up';

    btn_down = uibutton(fig, 'ButtonPushedFcn', @one_down);
    btn_down.Position = [0.02*w 0.45*h 0.07*w 0.05*h];
    btn_down.Text = 'One Down';
    
    ax = uiaxes(fig);
    ax.Position = [0.1*w 0.2*h 0.8*w 0.75*h];
    plot(ax, Wavenumber, Norm, 'b', 'LineWidth', 2)
    
    title(ax, title_text, 'interpreter', 'latex', 'FontSize', 18)
    xlabel(ax, x_text, 'interpreter', 'latex', 'FontSize', 14)
    ylabel(ax, y_text, 'interpreter', 'latex', 'FontSize', 14)
    
    ax.XLim = [Wavenumber(end) Wavenumber(1)];
    ax.YLim = [0 1];
    
    function sliderMoving(event,~)
        index = floor(event.Value);
        i = mod(index, x_size);
        if i == 0
            i = x_size;
        end
        j = ceil(index/x_size);
        
        Index = num2str(index); Index = strcat('Index= ', Index);
        I = num2str(i); 
        J = num2str(j);
        IJ = strcat('I= ', I, '; J= ', J);

        txt.Value = {Index, IJ};
    
        Wavenumber = sorted_data{i,j,1};
        Intensity = sorted_data{i,j,2};
        Norm = normalize(Intensity,1,'range');
        plot(ax, Wavenumber, Norm, 'b', 'LineWidth', 2)
    end

    function one_up(src,~)
        index = sscanf(txt.Value{1}, 'Index=%d');
        if index < ij
            index = index + 1;
        end

        sld.Value = index;

        i = mod(index, x_size);
        if i == 0
            i = x_size;
        end
        j = ceil(index/x_size);
        
        Index = num2str(index); Index = strcat('Index= ', Index);
        I = num2str(i); 
        J = num2str(j);
        IJ = strcat('I= ', I, '; J= ', J);

        txt.Value = {Index, IJ};
    
        Wavenumber = sorted_data{i,j,1};
        Intensity = sorted_data{i,j,2};
        Norm = normalize(Intensity,1,'range');
        plot(ax, Wavenumber, Norm, 'b', 'LineWidth', 2)
    end

    function one_down(src,~)
        index = sscanf(txt.Value{1}, 'Index=%d');
        if index > 1
            index = index - 1;
        end

        sld.Value = index;

        i = mod(index, x_size);
        if i == 0
            i = x_size;
        end
        j = ceil(index/x_size);
        
        Index = num2str(index); Index = strcat('Index= ', Index);
        I = num2str(i); 
        J = num2str(j);
        IJ = strcat('I= ', I, '; J= ', J);

        txt.Value = {Index, IJ};
    
        Wavenumber = sorted_data{i,j,1};
        Intensity = sorted_data{i,j,2};
        Norm = normalize(Intensity,1,'range');
        plot(ax, Wavenumber, Norm, 'b', 'LineWidth', 2)
    end
end