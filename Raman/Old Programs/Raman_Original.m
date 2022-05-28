data_name = get_data();
data = table2array(readtable(data_name));

image_name = get_image();
image = imread(image_name);

figure(), imshow(image)

[x_coord, y_coord] = ginput(2);

[i, j] = line_intersection(data, image, x_coord, y_coord);

figure()
for a = 1:length(i)
    [X, Y] = text_to_curve(data, i(a), j(a));
    hold on
    plot(X,Y)
end

figure()
[X, Y, Z] = curve_to_surface('line', data, i, j);
plot3(X, Z, Y)
xlabel('Wavelength (nm)')
ylabel('Transition')
zlabel('Intensity')

figure()
[X, Y, Z] = curve_to_surface('surface', data, i, j);
mesh(X, Z, Y)
xlabel('Wavelength (nm)')
ylabel('Transition')
zlabel('Intensity')