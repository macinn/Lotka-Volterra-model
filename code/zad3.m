data = readtable('HudsonBay.csv');
data = renamevars(data, data.Properties.VariableNames, ["t", "x", "y"]);
data.t = normalize(data.t, 'range');

N1 = 10;
N2 = 10;
N3 = 10;
N4 = 10;
J = zeros(N1,N2,N3, N4);
RX = linspace(0,10,N1);
RXY = linspace(-1,0,N2);
RXX = linspace(-2,0,N3);
X = linspace(min(data.x),max(data.x),N4);
dt = mean(diff(data.t));
T = 0 : dt : 1;
if true
for it1 = 1 : N1
    for it2 = 1 : N2
        for it3 = 1 : N3
            for it4 = 1 : N4
                rx = RX(it1);
                rxy = RXY(it2);
                rxx = RXX(it3);
                x = zeros(size(T));
                x(1) = X(it4);
                for i = 2 : length(T)
                    x1 = -(dt*rx - (dt^2*rx^2 + 2*dt^2*rx*rxy*data.y(i) + dt^2*rxy^2*data.y(i)^2 - 2*dt*rx - 2*dt*rxy*data.y(i) - 4*rxx*x(i-1)*dt + 1)^(1/2) + dt*rxy*data.y(i) - 1)/(2*dt*rxx);
                    x2 = -(dt*rx + (dt^2*rx^2 + 2*dt^2*rx*rxy*data.y(i) + dt^2*rxy^2*data.y(i)^2 - 2*dt*rx - 2*dt*rxy*data.y(i) - 4*rxx*x(i-1)*dt + 1)^(1/2) + dt*rxy*data.y(i) - 1)/(2*dt*rxx);
                    if isreal(x1) && x1>-10
                        x(i) = x1;
                    elseif isreal(x2) && x2>-10
                        x(i) = x2;
                    else
                        x(i) = Inf;
                        break
                    end
                end
            
                J(it1, it2, it3, it4) = sum((data.x'-x).*(data.x'-x));
            end
        end
    end
end
[~, ind] = min(J(:));
[it1, it2, it3, it4] = ind2sub(size(J), ind);
A = fminsearch(@(X) ApproxX(X(1), X(2), X(3), X(4)), [RX(it1), RXY(it2), RXX(it3), X(it4)]);
fprintf("r_x=%f; r_xy=%f; r_xx=%f; x0=%f\n", A(1), A(2), A(3), A(4));
r_x=A(1); r_xy=A(2); r_xx=A(3); x0=A(4);
end
J = zeros(N1,N2,N3, N4);
RY = linspace(0,10,N1);
RYX = linspace(0,1,N2);
RYY = linspace(-2,0,N3);
Y = linspace(min(data.y), max(data.y), N4);
for it1 = 1 : N1
    for it2 = 1 : N2
        for it3 = 1 : N3
            for it4 = 1 : N4
                ry = RY(it1);
                ryx = RYX(it2);
                ryy = RYY(it3);      
                y = zeros(size(T));
                y(1) = Y(it4);
                for i = 2 : length(T)
                    y1 = -(dt*ry - (dt^2*ry^2 + 2*dt^2*ry*ryx*data.x(i) + dt^2*ryx^2*data.x(i)^2 - 2*dt*ry - 2*dt*ryx*data.x(i) - 4*ryy*y(i-1)*dt + 1)^(1/2) + dt*ryx*data.x(i) - 1)/(2*dt*ryy);
                    y2 = -(dt*ry + (dt^2*ry^2 + 2*dt^2*ry*ryx*data.x(i) + dt^2*ryx^2*data.x(i)^2 - 2*dt*ry - 2*dt*ryx*data.x(i) - 4*ryy*y(i-1)*dt + 1)^(1/2) + dt*ryx*data.x(i) - 1)/(2*dt*ryy);
                    if isreal(y1) && y1>-10
                        y(i) = y1;
                    elseif isreal(y2) && y2>-10
                        y(i) = y2;
                    else
                        y(i) = Inf;
                        break
                    end
                end
                J(it1, it2, it3, it4) = sum((data.y'-y).*(data.y'-y));
            end
        end
    end
end
[v, ind] = min(J(:));
[it1, it2, it3, it4] = ind2sub(size(J), ind);
A = fminsearch(@(X) ApproxY(X(1), X(2), X(3), X(4)), [RY(it1), RYX(it2), RYY(it3), Y(it4)]);
fprintf("r_y=%f; r_yx=%f; r_yy=%f; y0=%f;\n", A(1), A(2), A(3), A(4));
r_y=A(1); r_yx=A(2); r_yy=A(3); y0=A(4);

J = @(X) fun3(X(1), X(2),X(3),X(4),X(5),X(6), X(7), X(8));
W = fminsearch(J, [x0, y0, r_x, r_y, r_xy, r_yx, r_xx, r_yy]);
disp(W)
hold on 
plot(data.t, data.x,'r');
plot(data.t, data.y,'b');
legend(["x" "y"])
f = @(t,x) [W(3)*x(1) + W(5)*x(1)*x(2)+ W(7)*x(1)*x(1); ...
    W(4)*x(2) + W(6)*x(1)*x(2) + W(8)*x(2)*x(2)];
[t,y] = ode45(f, [0 1], [W(1), W(2)]);
plot(t,y(:,1),'r--');
plot(t,y(:,2),'b--');

function [J] = ApproxX(rx ,rxy, rxx, x0)
data = readtable('HudsonBay.csv');
data = renamevars(data, data.Properties.VariableNames, ["t", "x", "y"]);
data(1:3,:)=[];
    data.t = normalize(data.t, 'range');
    dt = mean(diff(data.t));
    T = 0:dt:1;
    x = zeros(1, length(T));
    x(1) = x0;
    f = @(x,y) rx*x+rxy*x*y+rxx*x*x;
    for i = 2 : length(T)
        x1 = -(dt*rx - (dt^2*rx^2 + 2*dt^2*rx*rxy*data.y(i) + dt^2*rxy^2*data.y(i)^2 - 2*dt*rx - 2*dt*rxy*data.y(i) - 4*rxx*x(i-1)*dt + 1)^(1/2) + dt*rxy*data.y(i) - 1)/(2*dt*rxx);
        x2 = -(dt*rx + (dt^2*rx^2 + 2*dt^2*rx*rxy*data.y(i) + dt^2*rxy^2*data.y(i)^2 - 2*dt*rx - 2*dt*rxy*data.y(i) - 4*rxx*x(i-1)*dt + 1)^(1/2) + dt*rxy*data.y(i) - 1)/(2*dt*rxx);
        if isreal(x1) && x1>0
            x(i) = x1;
        elseif isreal(x2) && x2>0
            x(i) = x2;
        else
            x(i) = Inf;
            break
        end
    end

    J = sum((data.x'-x).*(data.x'-x));

end % function

function [J] = ApproxY(ry, ryx, ryy, y0)
data = readtable('HudsonBay.csv');
data = renamevars(data, data.Properties.VariableNames, ["t", "x", "y"]);
data(1:3,:)=[];
    data.t = normalize(data.t, 'range');   
    dt = mean(diff(data.t));
    T = 0:dt:1;
    y = zeros(1, length(T));
    y(1) = y0;
    for i = 2 : length(T)
        y1 = -(dt*ry - (dt^2*ry^2 + 2*dt^2*ry*ryx*data.x(i) + dt^2*ryx^2*data.x(i)^2 - 2*dt*ry - 2*dt*ryx*data.x(i) - 4*ryy*y(i-1)*dt + 1)^(1/2) + dt*ryx*data.x(i) - 1)/(2*dt*ryy);
        y2 = -(dt*ry + (dt^2*ry^2 + 2*dt^2*ry*ryx*data.x(i) + dt^2*ryx^2*data.x(i)^2 - 2*dt*ry - 2*dt*ryx*data.x(i) - 4*ryy*y(i-1)*dt + 1)^(1/2) + dt*ryx*data.x(i) - 1)/(2*dt*ryy);
        if isreal(y1) && y1>0
            y(i) = y1;
        elseif isreal(y2) && y2>0
            y(i) = y2;
        else
            y(i) = Inf;
            break
        end
    end
    J= sum((data.y'-y).*(data.y'-y));

end % function

function J = fun3(x1, y1, rx, ry, rxy, ryx, rxx, ryy)
data = readtable('HudsonBay.csv');
data = renamevars(data, data.Properties.VariableNames, ["t", "x", "y"]);
data(1:3,:)=[]; 
data.t = normalize(data.t, 'range');

    f = @(t, x) [rx*x(1) + rxy*x(1)*x(2)+rxx*x(1)*x(1); ...
        ry*x(2) + ryx*x(1)*x(2) + ryy*x(2)*x(2)];

    [t, y] = ode45(f, [0 1], [x1, y1]);
    X = interp1(t, y, data.t);
    J = sum( (X-[data.x, data.y]).*(X-[data.x, data.y]), 'all');
end % function