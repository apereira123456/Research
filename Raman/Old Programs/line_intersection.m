function [i, j] = line_intersection(data, image, x_coord, y_coord)
    [x_scans, y_scans, ~] = scan_size(data);
    [y_pixels, x_pixels, ~] = size(image);
    scale = [x_pixels/(x_scans-1), y_pixels/(y_scans-1)];
    
    x_pos = x_coord / scale(1);
    x_min = min(x_pos, [], 1);
    x_max = max(x_pos, [], 1);

    y_coord = y_pixels - y_coord;
    y_pos = y_coord / scale(2);
    y_min = min(y_pos, [], 1);
    y_max = max(y_pos, [], 1);

    x = floor(x_min) : 1 : ceil(x_max);
    y = floor(y_min) : 1 : ceil(y_max);

    m = (y_pos(2) - y_pos(1)) / (x_pos(2) - x_pos(1));
    b = (y_pos(1) - m * x_pos(1));

    hrz = @(y, m, b) [(y - b)./m; y];
    vrt = @(x, m, b) [x; m.*x + b];
        
    if m < 0
        h_int = ceil(hrz(y, m, b)');
        v_int = floor(vrt(x, m, b)');
        hv_int = [h_int; v_int];
    else
        h_int = ceil(hrz(y, m, b)');
        v_int = ceil(vrt(x, m, b)');
        hv_int = [h_int; v_int];
    end

    ex_bdry = find(hv_int(:,1) < floor(x_min) | hv_int(:,1) < 0);
    ex_bdry = [find(hv_int(:,1) > ceil(x_max) | hv_int(:,1) > x_scans - 1); ex_bdry];
    ex_bdry = [find(hv_int(:,2) < y_min | hv_int(:,2) < 0 ); ex_bdry];
    ex_bdry = [find(hv_int(:,2) > y_max | hv_int(:,2) > y_scans - 1); ex_bdry];
    hv_int(ex_bdry,:) = [];
        
    if m < 0
        unq = unique(hv_int,'rows');
        srtd = flip(sortrows(unq,[1 2],{'descend' 'ascend'}));
    else
        srtd = unique(hv_int,'rows');
    end

    i = srtd(:,1);
    j = srtd(:,2);
    a = [x_pos(1), x_pos(2)];
    b = [y_pos(1), y_pos(2)];

    figure()
    plot(a,b)
    xlim([0 x_scans-1])
    ylim([0 y_scans-1])
    set(gca,'xtick',0:1:x_scans-1)
    set(gca,'ytick',0:1:y_scans-1)
    grid on
end