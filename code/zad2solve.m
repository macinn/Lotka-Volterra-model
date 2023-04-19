W =[9.2205   -8.1276   -0.0977    0.0545    0.0001   -0.0098];
syms x y
 assume(x>0 & y>0)
eqn1 = 0 == W(1)*x + W(3)*x*y + W(5)*x*x;
eqn2 = 0 == W(2)*y + W(4)*x*y + W(6)*y*y;

[xx, yy] = solve([eqn1, eqn2]);
double(xx)
double(yy)
