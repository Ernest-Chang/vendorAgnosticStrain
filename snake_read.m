function [xs, ys] = snake_read(image)
%SNAKE_READ   Read snake control points
%   [XS, YS] = SNAKE_READ(IMAGE) displays the image in a new figure, then
%   records the position of each click of button 1 of the mouse in the
%   figure, and stops after recording the position of one click of another
%   button. The snake is drawn as it goes along. XS is a column vector of
%   x-coordinates and YS a column vector of y-coordinates.
%
%   See also SNAKE_DEMO

figure;
imagesc(image),colormap(gray);     % display image
hold on;           % and keep it there while we plot
xs = [];
ys = [];           % initialisations
xold = 0;
yold = 0;
but = 1;

while but == 1                      % while button 1 being pressed
    [xi, yi, but] = ginput(1);      % get a point
    xs = [xs; xi];                  % append its coords to vector
    ys = [ys; yi];
    if xold;
        plot([xold xi], [yold yi], 'go-');  % draw as we go
    else
        plot(xi, yi, 'go');         % first point on its own
    end;
    xold = xi;
    yold = yi;
end;

plot([xi xs(1)], [yi ys(1)], 'g-'); % join first to last points
hold off;

end
