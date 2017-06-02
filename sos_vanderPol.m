% % example for checking for checking lyapunov function in a region 

close all
clear 
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
Lagr1.degree = 4;
[Lagr1.poly, Lagr1.coeff] = polynomial(state.x,Lagr1.degree);

% Lagrangian2 for V_dot
Lagr2.degree = 4;
[Lagr2.poly, Lagr2.coeff] = polynomial(state.x,Lagr2.degree);

% decision variables 
% for finding lagrange polynomials 
var = [ Lyap.coeff; Lagr1.coeff;Lagr2.coeff];

r = 2; d = 0.01;
d_min = 1e-6;
d_max = 100;
% for i = 1:20
%     d = (d_min+d_max)/2;
    % Region 
    d = 3.9
    Reg.f = state.x(1)^2 + state.x(2)^2 - d;
    
    % constraints 
    F  = [sos( Lyap.poly + Lagr1.poly*Reg.f )
         sos( Lagr1.poly)
         sos(-Lyap.jaco_x*state.dx + Lagr2.poly*Reg.f)
         sos( Lagr2.poly)
         Lyap.coeff(1)==0
         ];

    [sol,v,Q] = solvesos(F,[],[],var); 
    flag = 1;
    if (sol.problem ~=0) % || (sum(abs(Q{1}(:))) - abs(Q{1}(1)) < 1e-3)
        flag = 0;
    end
    
    if flag == 1
        disp(['Feasible for d = ', num2str(d)]);
        d_min =d; 
    else
        disp(['Infeasible for d = ',num2str(d)])
        d_max = d;
    end
% end
if d_min < 1e-3
    disp(['Could not compute region of attraction']);
else
    disp(['Maximum feasible circle  with r^2  in  ',num2str((d_min)),' to ', num2str(d_max)]);
    disp(['radius is approximately',num2str(sqrt(d))])
end

%%
figure;
[x,y]= meshgrid(-5:0.1:5,-5:0.1:5);
Map = x*0;
V_grid = x*0;
dV_grid = x*0;
for i =1:length(x(:))
    V_grid(i) = value(replace(Lyap.poly,state.x,[x(i),y(i)]));
    dV_grid(i) = value(replace(Lyap.jaco_x*state.dx,state.x,[x(i),y(i)]));
    if((V_grid(i)>=0) && (dV_grid(i)<0))
        Map (i) =100;
    end
end
V_grid = reshape(V_grid,size(x));
dV_grid = reshape(dV_grid,size(x));
surf(x,y,V_grid);
figure;
surf(x,y,dV_grid);
figure;
surf(x,y,Map);
view(0,90);