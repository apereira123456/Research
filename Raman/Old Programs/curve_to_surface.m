function [X, Y, Z] = curve_to_surface(mode, data, i, j)
    [~, ~, s_size] = scan_size(data);
    
    n = size(i);
    [X, Y, Z] = deal(zeros(s_size(1), n(1)), zeros(s_size(1), n(1)), zeros(s_size(1), n(1)));
    
    if strcmp(mode,'line')
        for k = 1:n(1)
            [Wavenumber, Intensity] = text_to_curve(data, i(k), j(k));
            X(:,k) = Wavenumber;
            Y(:,k) = Intensity;
            Z(:,k) = k;
        end
        
    elseif strcmp(mode,'surface')
        for k = 1:n(1)
            [Wavenumber, Intensity] = text_to_curve(data, i(k), j(k));
            X(:,k) = Wavenumber;
            Y(:,k) = Intensity;
        end

        [X, Z] = meshgrid(X(:,1), 20*i(:,1));
        Y = (Y)';
    else
        disp('Incompatible mode')
    end
end