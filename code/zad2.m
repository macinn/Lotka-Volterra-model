data = readtable('dane21.csv');
J = @(X) fun(X(1), X(2),X(3),X(4),X(5),X(6));

r_x=6.136346; r_xy=-0.069893; r_xx=-0.000000;
r_y=-5.884900; r_yx=0.059108; r_yy=-0.033141;

[W, val] = fminsearch(J, [r_x, r_y, r_xy, r_yx, r_xx, r_yy]);
disp(W);

f = figure;
hold on 
f.Position = [100 100 1000 350];
plot(data.t, data.x,'r');
plot(data.t, data.y,'b');
f = @(t,x) [W(1)*x(1) + W(3)*x(1)*x(2)+ W(5)*x(1)*x(1); ...
    W(2)*x(2) + W(4)*x(1)*x(2) + W(6)*x(2)*x(2)];
[t,y] = ode45(f, [0 3], [data.x(1), data.y(1)]);
plot(t,y(:,1),'r--');
plot(t,y(:,2),'b--');

function J = fun(rx, ry, rxy, ryx, rxx, ryy)
    data = readtable('dane21.csv');
    
    f = @(t, x) [rx*x(1) + rxy*x(1)*x(2)+rxx*x(1)*x(1); ...
        ry*x(2) + ryx*x(1)*x(2) + ryy*x(2)*x(2)];

    [t, y] = ode45(f, [0 3], [data.x(1), data.y(1)]);
    X = interp1(t, y, data.t);
    J = sum( (X-[data.x, data.y]).*(X-[data.x, data.y]), 'all');
end