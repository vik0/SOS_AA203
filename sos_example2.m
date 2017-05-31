% Example for checking for lyapunov function for a simple linear system

close all
clear 
clc

% lets try a simple system of form x_dot = -x
x = sdpvar(2,1);
[p,c] = polynomial(x,2);
j = jacobian(p,x);

% dynamics of system 
% x_dot = A*x;

A = [-1, 0;
      0,-1];
x_dot = A*x;
F = [sos(p), sos(-j*x_dot)];

[sol,v,Q] = solvesos(F,[],[],c)
disp('Monomials used are ');
sdisplay(v{1}')
disp('PSD Matrix is ');
Q{1}
disp('we can see that a PSD matrix (not constant) EXISTS');

disp('Try changing the matrix A and see that PSD matrix will be');
disp('constant in certain cases telling that lyapunov function doesnot exist');