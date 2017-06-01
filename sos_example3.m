% example for checking for checking lyapunov function in a region 

close all
clear all
clc

% state dynamics
A = -eye(2);
dyan_simple = @(x)A*x;

% Vanderpol dynamics (inverted in time)
% dot(x) = -y
% dot(y) = -(mu(1-x^2)y - x)
mu = 0.01;
dyan_vp = @(x) [- x(2);
    -mu*(1-x(1)^2)*x(2) + x(1)];

% state
state.dim = 2;
state.x = sdpvar(state.dim,1);
% state.dx = dyan_simple(state.x);
state.dx = dyan_vp(state.x);


% Lyapunov - p
Lyap.degree = 6;
[Lyap.poly, Lyap.coeff] = polynomial(state.x,Lyap.degree); 
Lyap.jaco_x = jacobian(Lyap.poly,state.x); 

% Region 
d = sdpvar(1,1);
Reg.poly = sum(state.x.^2) - d; 

% Lagrangian1
Lagr1.degree = 4;
[Lagr1.poly, Lagr1.coeff] = polynomial(state.x,Lagr1.degree);

% Lagrangian2 for V_dot
Lagr2.degree = 4;
[Lagr2.poly, Lagr2.coeff] = polynomial(state.x,Lagr2.degree);


% constraints 
F  = [sos( Lyap.poly + Lagr1.poly*Reg.poly )
     sos( Lagr1.poly)
     sos(-Lyap.jaco_x*state.dx + Lagr2.poly*Reg.poly)
     sos( Lagr2.poly)];

% decision variables 
% for finding lagrange polynomials 
var = [ Lyap.coeff; Lagr1.coeff;Lagr2.coeff];

d_min = 1e-6;
d_max = 100;
for i = 1:20
    d = (d_min+d_max)/2;
    [sol,v,Q] = solvesos(F,[],[],var);
    flag = 1;
    if (sol.problem ~=0) || (sum(abs(Q{1 }(:))) - abs(Q{1}(1)) < 1e-6)...
            || (sum(abs(Q{3}(:))) < 1e-6)
        flag = 0;
    end
    
    if flag == 1
        disp(['Feasible for d = ', num2str(d)]);
        d_min =d;
    else
        disp(['Infeasible for d = ',num2str(d)])
        d_max = d;
    end
end
if d_min < 1e-3
    disp(['Could not compute region of attraction']);
else
    disp(['Maximum circle  with r^2  in  ',num2str((d_min)),' to ', num2str(d_max)]);
end