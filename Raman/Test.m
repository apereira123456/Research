x = linspace(900, 1900);
y1 = randn(size(x)) + cumsum(ones(size(x)))*0.3 + 1;
y2 = randn(size(x)) + cumsum(ones(size(x)))*0.9 + 2;
y1d = detrend(y1);
y2d = detrend(y2);
figure(1)
subplot(2,1,1)
plot(x, y1,    x, y2)
title('Original')
grid
subplot(2,1,2)
plot(x, y1d,    x, y2d)
title('Detrended')
grid