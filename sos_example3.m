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
mu = 10;
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
Lagr1.degree = 6;
[Lagr1.poly, Lagr1.coeff] = polynomial(state.x,Lagr1.degree);

% Lagrangian2 for V_dot
Lagr2.degree = 6;
[Lagr2.poly, Lagr2.coeff] = polynomial(state.x,Lagr2.degree);

% decision variables 
% for finding lagrange polynomials 
var = [ Lyap.coeff; Lagr1.coeff;Lagr2.coeff];

options = sdpsettings('sos.newton',1,'sos.congruence',1,'verbose',0);

r = 2; d = 0.01;
d_min = 1e-6;
d_max = 10;
for i = 1:10
    d = (d_min+d_max)/2;
    % Region 
    Reg.f = state.x(1)^2 + state.x(2)^2 - d;
    
    % constraints 
    F  = [sos( Lyap.poly + Lagr1.poly*Reg.f )
         sos( Lagr1.poly)
         sos(-Lyap.jaco_x*state.dx + Lagr2.poly*Reg.f)
         sos( Lagr2.poly)
         Lyap.coeff(1)==0
         sum(Lyap.coeff) == 1
         ];

    [sol,v,Q] = solvesos(F,[],options,var); 
    flag = 1;
    if (sol.problem ~=0) || (sum(abs(Q{1}(:))) - abs(Q{1}(1)) < 1e-3)
        flag = 0;
    end
    
    if flag == 1
        disp(['Feasible for d = ', num2str(d)]);
        V_c0 = value(Lyap.coeff);
        d_min =d; 
    else
        disp(['Infeasible for d = ',num2str(d)])
        d_max = d;
    end
end
if d_min < 1e-3
    disp(['Could not compute region of attraction']);
else
    disp(['Maximum feasible circle  with r^2  in  ',num2str((d_min)),' to ', num2str(d_max)]);
    disp(['radius is approximately',num2str(sqrt(d))])
end

%% 
delta = 0.1;
[x1,y1]= meshgrid(-5:delta:5,-5:delta:5);
for i = 1:length(x1(:))
    X(i,:) = [x1(i), y1(i)];
    dX(i,:) = dyan_vp(X(i,:));
end
V = V_c0(2)*X(:,1)+V_c0(3)*X(:,2)+V_c0(4)*X(:,1).^2+V_c0(5)*X(:,1).*X(:,2)+V_c0(6)*X(:,2).^2 ... ;
    +X(:,1).^3*V_c0(7)+X(:,1).^2 .*X(:,2)*V_c0(8)+X(:,1).*X(:,2).^2*V_c0(9)+X(:,2).^3*V_c0(10)+ ...
    X(:,1).^4*V_c0(11)+X(:,1).^3 .*X(:,2)*V_c0(12)+X(:,1).^2 .*X(:,2).^2*V_c0(13)+ ...
    X(:,1).*X(:,2).^3*V_c0(14)+X(:,2).^4*V_c0(15);

dVdX(:,1) = V_c0(2) + 2*V_c0(4)*X(:,1) + V_c0(5)*X(:,2)  ...;
       +3*X(:,1).^2*V_c0(7) + 2*X(:,1).*X(:,2)*V_c0(8) + X(:,2).^2*V_c0(9)+ ...
       4*X(:,1).^3*V_c0(11) + 3*X(:,1).^2 .*X(:,2)*V_c0(12) + 2*X(:,1).*X(:,2).^2*V_c0(13)+ ...
       X(:,2).^3*V_c0(14);
dVdX(:,2) = V_c0(3) + V_c0(5)*X(:,1) + 2*V_c0(6)*X(:,2)  ...;
       + X(:,1).^2 *V_c0(8)+ 2* X(:,1).*X(:,2)*V_c0(9)+3*X(:,2).^2*V_c0(10)+ ...
       X(:,1).^3 *V_c0(12)+2*X(:,1).^2 .*X(:,2)*V_c0(13)+ ...
       3*X(:,1).*X(:,2).^2*V_c0(14)+4*X(:,2).^3*V_c0(15);

dV = sum(dVdX.*dX,2);
%%
Map = x1*0+1;
Map(dV>=0) = 0; 
Map(V<=0)=0; 
V = reshape(V,size(x1));
dV = reshape(dV,size(x1));

dyan_tp = @(t, x) [ x(2);
mu*(1-x(1)^2)*x(2) - x(1)];
[t, a] = ode45(dyan_tp, [1, 20], [0, 1]);

figure;
% surf(x1,y1,V,'EdgeColor','none');
surf(x1,y1,V);
title('Lyapunov');
axis([-5 5 -5 5 0 3.9]);
hold on;
plot(a(:,1),a(:,2),'r','LineWidth',2);
view(0,90);

figure;
surf(x1,y1,dV);
title('Grad Lyapunov');
axis([-5 5 -5 5 -10 0])


figure;
surf(x1,y1,Map); hold on;
title('Region with V>0 and dV<0');
plot(a(:,1),a(:,2));
view(0,-90);

figure(1);
