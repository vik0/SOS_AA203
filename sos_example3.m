% example for checking for checking lyapunov function in a region 

close all
clear all
clc

% lets try a simple system of form x_dot = -x
x = sdpvar(2,1);
[p,c] = polynomial(x,2);
j = jacobian(p,x);

% dynamics of system 
% x_dot = 
% for vanderpol oscilltor (inverted in time)
% dot(x) = -y
% dot(y) = -(mu(1-x^2)y - x)

mu = 0.01;
x_dot =     [ -x(2);
    -mu*(1-x(1)^2)*x(2)+x(1)];


d = sdpvar(1,1);
reg = d^2 - (x(1)^2 + x(2)^2) ;
% F = [sos(p), sos(-j*x_dot),sos(reg)];
F = [sos(p), sos(-j*x_dot)];

% [sol,v,Q] = solvesos(F,d,[],[c;d])
[sol,v,Q] = solvesos(F,[],[],[c])
disp('Monomials used are ');
sdisplay(v{1}')
disp('PSD Matrix is ');
Q{1}

disp('if  PSD matrix is  constant then lyapunov function doesnot exist');