data = readtable('dane21.csv');
N1 = 20; N2 = 20; N3 = 20;
J = zeros(N1,N2,N3);
RX = linspace(-100,100,N2);
RXY = linspace(-1,1,N2);
RXX = linspace(-1,1,N3);
dt = mean(diff(data.t));
T = 0 :  dt : 3;

for it1 = 1 : N1
    for it2 = 1 : N2
        for it3 = 1 : N3
            rx = RX(it1);
            rxy = RXY(it2);
            rxx = RXX(it3);
            x = zeros(size(T));
            x(1) = data.x(1);
            f = @(x,y) rx*x+rxy*x*y+rxx*x*x;
%             Metoda Eulera
%             for i = 2 : length(T)
%                 x(i) = x(i-1) + dt * f(x(i-1), data.y(i-1));
%             end

%             Adams-Bashforth
            x(2) = x(1) + dt * f(x(1), data.y(1));
            x(3) = x(2) + dt * f(x(2), data.y(2));
            for i = 4 : length(T)
                x(i) = x(i-1) + dt/12 * (23*f(x(i-1), data.y(i-1))...
                    -16*f(x(i-2), data.y(i-2)) +5*f(x(i-3), data.y(i-3)));
            end

%             Niejawna metoda Eulera
%             for i = 2 : length(T)
%                 x1 = -(dt*rx - (dt^2*rx^2 + 2*dt^2*rx*rxy*data.y(i) + dt^2*rxy^2*data.y(i)^2 - 2*dt*rx - 2*dt*rxy*data.y(i) - 4*rxx*x(i-1)*dt + 1)^(1/2) + dt*rxy*data.y(i) - 1)/(2*dt*rxx);
%                 x2 = -(dt*rx + (dt^2*rx^2 + 2*dt^2*rx*rxy*data.y(i) + dt^2*rxy^2*data.y(i)^2 - 2*dt*rx - 2*dt*rxy*data.y(i) - 4*rxx*x(i-1)*dt + 1)^(1/2) + dt*rxy*data.y(i) - 1)/(2*dt*rxx);
%                 if isreal(x1) && isreal(x2) && (x1>0 || x2>0)
%                     x(i) = max(x1, x2);
%                 else
%                     x(i) = Inf;
%                     break
%                 end
%             end

            J(it1, it2, it3) = sum((data.x'-x).*(data.x'-x));
        end
    end
end
[v, ind] = min(J(:));
[it1, it2, it3] = ind2sub(size(J), ind);
A = fminsearch(@(X) ApproxX(X(1), X(2), X(3)), [RX(it1), RXY(it2), RXX(it3)]);
fprintf("r_x=%f; r_xy=%f; r_xx=%f;\n", A(1), A(2), A(3));

f = @(x,y) A(1)*x+A(2)*x*y+A(3)*x*x;
x(2) = x(1) + dt * f(x(1), data.y(1));
x(3) = x(2) + dt * f(x(2), data.y(2));
for i = 4 : length(T)
    x(i) = x(i-1) + dt/12 * (23*f(x(i-1), data.y(i-1))...
        -16*f(x(i-2), data.y(i-2)) +5*f(x(i-3), data.y(i-3)));
end
f=figure;
hold on
plot(data.t, data.x, 'r');
plot(T, x, 'r--');
f.Position=[100 100 1000 350];

function [J] = ApproxX(rx ,rxy, rxx)
data = readtable('dane21.csv');
dt = mean(diff(data.t));
T = 0 :  dt : 3;
x = zeros(1, length(T));
x(1) = data.x(1);

for i = 2 : length(T)
    x1 = -(dt*rx - (dt^2*rx^2 + 2*dt^2*rx*rxy*data.y(i) + dt^2*rxy^2*data.y(i)^2 - 2*dt*rx - 2*dt*rxy*data.y(i) - 4*rxx*x(i-1)*dt + 1)^(1/2) + dt*rxy*data.y(i) - 1)/(2*dt*rxx);
    x2 = -(dt*rx + (dt^2*rx^2 + 2*dt^2*rx*rxy*data.y(i) + dt^2*rxy^2*data.y(i)^2 - 2*dt*rx - 2*dt*rxy*data.y(i) - 4*rxx*x(i-1)*dt + 1)^(1/2) + dt*rxy*data.y(i) - 1)/(2*dt*rxx);
    if isreal(x1) && isreal(x2) && (x1>0 || x2>0)
        x(i) = max(x1, x2);
    else
        x(i) = Inf;
        break
    end
end
J = sum((data.x'-x(1,:)).*(data.x'-x(1, :)));

end % function
