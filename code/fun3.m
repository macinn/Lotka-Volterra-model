function J = fun3(x1, y1, rx, ry, rxy, ryx, rxx, ryy)
data = readtable('IsleRoyale.csv');
data = renamevars(data, data.Properties.VariableNames, ["t", "x", "y"]);
data.t = normalize(data.t, 'range');

    f = @(t, x) [rx*x(1) + rxy*x(1)*x(2)+rxx*x(1)*x(1); ...
        ry*x(2) + ryx*x(1)*x(2) + ryy*x(2)*x(2)];

    [t, y] = ode45(f, [0 1], [x1, y1]);
    X = interp1(t, y, data.t);
    J = sum( (X-[data.x, data.y]).*(X-[data.x, data.y]), 'all');
end