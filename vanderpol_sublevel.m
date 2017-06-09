% SUBLEVEL SETS 

close all
clear
clc

mu = 1;
vander = @(x) [-x(2); -mu*(1-x(1)^2)*x(2) + x(1)];

x = sdpvar(2,1);
dx = vander(x);

% Lyap
[V_p,V_c] = polynomial(x,4);
dV_p = jacobian(V_p,x);

% Lagrangian 
[l_p,l_c] = polynomial(x,4);

% for initialization 
V_c0 = zeros(15,1);
V_c0(1:6)  = 0.5*[0; 0; 0; 1.5; -0.5 ; 1];
V_p0  = replace(V_p, V_c,V_c0);
dV_p0 = replace(dV_p,V_c,V_c0);
r0 = 0.5;

max_iter2 = 100;
for i = 1:max_iter2
    % binary search 
    r_max = 10;
    r_min = r0;
    l_c0 = [];
    iter = 1; max_iter = 10;
    while(iter<max_iter)
        r0 = (r_min +r_max)/2;
        % constraints
        F = [sos(l_p), sos(sum(x.^2)*1e-6-dV_p0*dx + l_p*(V_p0-r0))];
        var = l_c;
        options = sdpsettings('sos.newton',1,'sos.congruence',1,'verbose',0);
        [sol,v,Q] = solvesos(F,[],options,var);
        if sol.problem ~=0
            r_max = r0;
        else
            r_min = r0;
            l_c0 = value(l_c);
        end
         iter = iter + 1;
    end
    disp(r0)
    if length(l_c0)==0
        break;
    end
    
    % find Lyap
    sdpvar r
    l_p0 = replace(l_p,l_c,l_c0);
    F2 = [
        % sos(V_p)
        sos( - dV_p*dx + l_p0*(V_p-r))
        % r>=0  
        % V_c(3)+V_c(5)==1
        sum(V_c)==1
        V_c(1)==0
        ];
    var2 = [V_c; r];
    [sol2,v2,Q2] = solvesos(F2,-r,options,var2);
    
    if sol2.problem == 0
        V_c0 = value(V_c);
        V_p0  = replace(V_p, V_c,V_c0);
        dV_p0 = replace(dV_p,V_c,V_c0);
        r0 = value(r);
    else
        break;
    end
    disp(r0);
end

%% 
delta = 0.1;
[x1,y1]= meshgrid(-5:delta:5,-5:delta:5);
for i = 1:length(x1(:))
    X(i,:) = [x1(i), y1(i)];
    dX(i,:) = vander(X(i,:));
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

Map = x1*0+1; Map(V<=0)=0;Map(dV>=0) = 0;

V = reshape(V,size(x1));
dV = reshape(dV,size(x1));

dyan_tp = @(t, x) [ x(2);
mu*(1-x(1)^2)*x(2) - x(1)];
[t, a] = ode45(dyan_tp, [1, 20], [0, 1]);
figure;
surf(x1,y1,V);
% title('Lyapunov');
axis([-5 5 -5 5 0 r0]); hold on;
plot(a(:,1),a(:,2),'r','LineWidth',2);
view(0,90);

% 
% figure;
% surf(x1,y1,dV);
% title('Grad Lyapunov');
% axis([-5 5 -5 5 -10 0])
%  
% figure;;
% surf(x1,y1,Map);
% title('Region with V>0 and dV<0');