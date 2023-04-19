function J = fun(rx, ry, rxy, ryx, rxx, ryy)
    data = readtable('dane21.csv');
    
    f = @(t, x) [rx*x(1) + rxy*x(1)*x(2)+rxx*x(1)*x(1); ...
        ry*x(2) + ryx*x(1)*x(2) + ryy*x(2)*x(2)];

    [t, y] = ode45(f, [0 3], [data.x(1), data.y(1)]);
    X = interp1(t, y, data.t);
    J = sum( (X-[data.x, data.y]).*(X-[data.x, data.y]), 'all');
end