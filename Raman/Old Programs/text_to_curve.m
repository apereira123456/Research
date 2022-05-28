function [Wavenumber, Intensity] = text_to_curve(data, i, j)   
    [x_size, ~, s_size] = scan_size(data);
    a = 1;
    Wavenumber = zeros(s_size,1);
    Intensity = zeros(s_size,1);

    while a <= s_size(1)
        n = a + i*s_size(1) + j*x_size(1)*s_size(1);
        Wavenumber(a,1) = data(a,3);
        Intensity(a,1) = data(n,4);
        a = a + 1;
    end
end