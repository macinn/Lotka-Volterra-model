    data = readtable('dane21.csv');
    N1 = 50;
    N2 = 50;
    N3 = 50;
    J = zeros(N1,N2,N3);
    RY = linspace(-100,100,N1);
    RYX = linspace(-1,1,N2);
    RYY = linspace(-1,1,N3);
    dt = mean(diff(data.t));
    T = 0 :  dt : 3;
    for it1 = 1 : N1
        for it2 = 1 : N2
            for it3 = 1 : N3
                ry = RY(it1);
                ryx = RYX(it2);
                ryy = RYY(it3);
                y = zeros(size(T));
                y(1) = data.y(1);
                f = @(x,y) ry*y+ryx*x*y+ryy*y*y;
    %             Metoda Eulera
    %             for i = 2 : length(T)
    %                 y(i) = y(i-1) + dt*f(data.x(i-1), y(i-1));
    %             end
    
    %             Metoda Adamsa-Bashfortha
                y(2) = y(1) + dt * f(data.x(1), y(1));
                y(3) = y(2) + dt * f(data.x(2), y(2));
                for i = 4 : length(T)
                    y(i) = y(i-1) + dt/12 * (23*f(data.x(i-1), y(i-1))...
                        -16*f(data.x(i-2), y(i-2)) +5*f(data.x(i-3), y(i-3)));
                end
    
    %             Niejawna metoda Eulera
    %             for i = 2 : length(T)
    %                 y1 = -(dt*ry - (dt^2*ry^2 + 2*dt^2*ry*ryx*data.x(i) + dt^2*ryx^2*data.x(i)^2 - 2*dt*ry - 2*dt*ryx*data.x(i) - 4*ryy*y(i-1)*dt + 1)^(1/2) + dt*ryx*data.x(i) - 1)/(2*dt*ryy);
    %                 y2 = -(dt*ry + (dt^2*ry^2 + 2*dt^2*ry*ryx*data.x(i) + dt^2*ryx^2*data.x(i)^2 - 2*dt*ry - 2*dt*ryx*data.x(i) - 4*ryy*y(i-1)*dt + 1)^(1/2) + dt*ryx*data.x(i) - 1)/(2*dt*ryy);
    %                 if isreal(y1) && isreal(y2) && (y1>0 || y2>0)
    %                     y(i) = max(y1, y2);
    %                 else
    %                     y(i) = Inf;
    %                     break
    %                 end
    %             end
    
                J(it1, it2, it3) = sum((data.y'-y).*(data.y'-y));
            end
        end
    end
    
    [v, ind] = min(J(:));
    [it1, it2, it3] = ind2sub(size(J), ind);
    A = fminsearch(@(X) ApproxY(X(1), X(2), X(3)), [RY(it1), RYX(it2), RYY(it3)]);
    fprintf("r_y=%f; r_yx=%f; r_yy=%f;\n", A(1), A(2), A(3));
    
    f = @(x,y) A(1)*y+A(2)*x*y+A(3)*y*y;
    y(2) = y(1) + dt * f(data.x(1), y(1));
    y(3) = y(2) + dt * f(data.x(2), y(2));
    for i = 4 : length(T)
        y(i) = y(i-1) + dt/12 * (23*f(data.x(i-1), y(i-1))...
            -16*f(data.x(i-2), y(i-2)) +5*f(data.x(i-3), y(i-3)));
    end
    
    f=figure;
    hold on
    plot(data.t, data.y, 'b');
    plot(T, y, 'b--');
    f.Position=[100 100 1000 350];
    
    function [J] = ApproxY(ry, ryx, ryy)
    data = readtable('dane21.csv');
    dt = mean(diff(data.t));
    T = 0 :  dt : 3;
    y = zeros(size(T));
    y(1) = data.y(1);
    for i = 2 : length(T)
        y1 = -(dt*ry - (dt^2*ry^2 + 2*dt^2*ry*ryx*data.x(i) + dt^2*ryx^2*data.x(i)^2 - 2*dt*ry - 2*dt*ryx*data.x(i) - 4*ryy*y(i-1)*dt + 1)^(1/2) + dt*ryx*data.x(i) - 1)/(2*dt*ryy);
        y2 = -(dt*ry + (dt^2*ry^2 + 2*dt^2*ry*ryx*data.x(i) + dt^2*ryx^2*data.x(i)^2 - 2*dt*ry - 2*dt*ryx*data.x(i) - 4*ryy*y(i-1)*dt + 1)^(1/2) + dt*ryx*data.x(i) - 1)/(2*dt*ryy);
        if isreal(y1) && isreal(y2) && (y1>0 || y2>0)
            y(i) = max(y1, y2);
        else
            y(i) = Inf;
            break
        end
    end
    J= sum((data.y'-y).*(data.y'-y));
    
    end % function