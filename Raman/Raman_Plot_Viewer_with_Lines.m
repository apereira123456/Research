function Raman_Plot_Viewer_with_Lines
    %% Clear
    close all; clear all; clc;

    % Specify interval over which to search for peak maximum
    sample_material = 'B4C';
    % sample_material = 'SiC';
        
    if strcmp(sample_material, 'B4C')
        peak_id = {'480', '530', '730', '1080'};
        search_loc = [460,500; 510,550; 700,760; 1040,1120];
    elseif strcmp(sample_material, 'SiC')
        peak_id = {'790', '970'};
        search_loc = [770,810; 950,990];
    end

    %% Data Selection Prompt
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
    sorted_data(:,:,1) = mat2cell(data(s_size:-1:1,3),s_size,1); 
    
    for j = 1:y_size
        for i = 1:x_size
            l = (i-1)*s_size + (j-1)*x_size*s_size + 1;
            b = i*s_size + (j-1)*x_size*s_size;

            sorted_data(i,j,2) = mat2cell(normalize(data(b:-1:l,4),1,'range'),s_size,1);
        end
    end
    
    %% Graph Labels
    title_text = '\textbf{Raman Viewer}';
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
    txt.Position = [0.9*w 0.5*h 0.07*w 0.05*h];
    txt.Value = {'Index=1', 'I=1; J=1'};

    % Suppress warning for fixed slider height
    id = 'MATLAB:ui:Slider:fixedHeight';
    warning('off',id)

    sld = uislider(fig, 'ValueChangingFcn', @(sld,event) sliderMoving(event,sld));
    sld.Position = [0.1*w 0.1*h 0.8*w 0.05*h];
    sld.Value = 1;
    sld.Limits = [1 ij];
    
    btn_up = uibutton(fig, 'ButtonPushedFcn', @one_up);
    btn_up.Position = [0.01*w 0.55*h 0.08*w 0.05*h];
    btn_up.Text = 'One Up';

    btn_down = uibutton(fig, 'ButtonPushedFcn', @one_down);
    btn_down.Position = [0.01*w 0.45*h 0.08*w 0.05*h];
    btn_down.Text = 'One Down';

    btn_save = uibutton(fig, 'ButtonPushedFcn', @save_fig);
    btn_save.Position = [0.01*w 0.1*h 0.08*w 0.05*h];
    btn_save.Text = 'Save Figure';
    
    ax = uiaxes(fig);
    ax.Position = [0.1*w 0.15*h 0.8*w 0.8*h];

    hold(ax,'on')

    for l = 1:length(peak_id)
        plot(ax, [search_loc(l,1) search_loc(l,1)], [0 1], 'r', 'LineWidth', 2)
        plot(ax, [search_loc(l,2) search_loc(l,2)], [0 1], 'r', 'LineWidth', 2)
    end

    pt = plot(ax, Wavenumber, Intensity, 'b', 'LineWidth', 2);
    
    title(ax, title_text, 'interpreter', 'latex', 'FontSize', 18)
    xlabel(ax, x_text, 'interpreter', 'latex', 'FontSize', 14)
    ylabel(ax, y_text, 'interpreter', 'latex', 'FontSize', 14)
    
    ax.XLim = [Wavenumber(1) Wavenumber(end)];
    ax.YLim = [0 1];
    
    function sliderMoving(event,~)
        delete(pt)

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
        pt = plot(ax, Wavenumber, Intensity, 'b', 'LineWidth', 2);
    end

    function one_up(src,~)
        delete(pt)

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
        pt = plot(ax, Wavenumber, Intensity, 'b', 'LineWidth', 2);
    end

    function one_down(src,~)
        delete(pt)

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
        pt = plot(ax, Wavenumber, Intensity, 'b', 'LineWidth', 2);
    end

    function save_fig(src,~)
        exportgraphics(ax,fullfile(path,'Raman_Plot.png'),'Resolution',300)
    end
end