% % example for checking for checking lyapunov function in a region
% ITERATIVELY

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
Lyap.degree = 4;
[Lyap.poly, Lyap.coeff] = polynomial(state.x,Lyap.degree); 
Lyap.jaco_x = jacobian(Lyap.poly,state.x); 

% Lagrangian1
Lagr1.degree = 2;
[Lagr1.poly, Lagr1.coeff] = polynomial(state.x,Lagr1.degree);

% Lagrangian2 for V_dot
Lagr2.degree = 2;
[Lagr2.poly, Lagr2.coeff] = polynomial(state.x,Lagr2.degree);


d = 1;
feasible_d = 0;
for i = 1:20
    % Region 
    Reg.f = state.x(1)^2 + state.x(2)^2 - d;
    
    % constraints 
    F1  = [sos( Lyap.poly + Lagr1.poly*Reg.f )
         sos( Lagr1.poly)
         sos(-Lyap.jaco_x*state.dx + Lagr2.poly*Reg.f)
         sos( Lagr2.poly)];
         
    % decision variables 
    % for finding lagrange polynomials 
    var1 = [ Lyap.coeff; Lagr1.coeff;Lagr2.coeff];

    [sol1,v,Q1] = solvesos(F1,[],[],var1); 
    if (sol1.problem ~=0) || (sum(abs(Q1{1}(:))) - abs(Q1{1}(1)) < 1e-3)
        break;
    end
    
    sdpvar t
    % Region 
    Reg.f = state.x(1)^2 + state.x(2)^2 - t;
    
    % constraints 
    F2  = [sos( Lyap.poly + Lagr1.poly*Reg.f ), sos(-Lyap.jaco_x*state.dx + Lagr2.poly*Reg.f), t<=100, t>=0];
         
    % decision variables 
    % for finding lagrange polynomials 
    var2 = [ Lyap.coeff; t];
    [sol2,v,Q2] = solvesos(F2,-t,[],var2); 
    if (sol2.problem ~=0) || (sum(abs(Q2{1}(:))) - abs(Q2{1}(1)) < 1e-3)
        break;
    end
    feasible_d = value(t);
    d = value(t);
    
end
disp(['Maximum circle  with r^2  is  ', num2str(feasible_d)]);