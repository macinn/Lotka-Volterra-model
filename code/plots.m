data = readtable('dane21.csv');


hold on 
plot(data.t, data.x);
plot(data.t, data.y);
f = @(t,x) [W(1)*x(1) + W(3)*x(1)*x(2)+ W(5)*x(1)*x(1); ...
    W(2)*x(2) + W(4)*x(1)*x(2) + W(6)*x(2)*x(2)];
[t,y] = ode45(f, [0 3], [data.x(1), data.y(1)]);
plot(t,y);