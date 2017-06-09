close all
clear 
clc

sdpvar x1 y1 x2 y2
x = [x1; y1; x2; y2];

[V_p,V_c] = polynomial(x,2,1);
dV_p = jacobian(V_p,x);

r0 = 0.5;
reg = (x1-x2)^2 + (y1-y2)^2 - r0^2;

[u_p1,u_c1] = polynomial(x,1);
[u_p2,u_c2] = polynomial(x,1);

[l1_p,l1_c] = polynomial(x,2);
[l2_p,l2_c] = polynomial(x,2);
[l3_p,l3_c] = polynomial(x,2);
[l4_p,l4_c] = polynomial(x,2);

for i=1:10
    u_c10 = randn(size(u_c1));
    u_c20 = randn(size(u_c2));

    u_p10 = replace(u_p1,u_c1,u_c10);
    u_p20 = replace(u_p2,u_c2,u_c20);
    dx0 = [u_p10; u_p20; x1-x2; y1-y2];
    F1 = [
        sos(V_p)
        %- l1_p*reg)
        %sos(l1_p)
        %sos(-V_p + l2_p*reg)
        sos(l2_p)
        sos(-dV_p*dx0 -l3_p*reg)
        sos(l3_p)
        sos(dV_p*dx0 +l4_p*reg)
        sos(l4_p)
        ];

    var1 = [V_c;l1_c;l2_c;l3_c;l4_c];

    options = sdpsettings('sos.newton',1,'sos.congruence',1,'verbose',0);

    [sol1,v,Q] = solvesos(F1,sum(V_c),options,var1);
    sol1

    V_c0 = value(V_c);
    dV_p0 = replace(dV_p,V_c,value(V_c));
    l3_p0 = replace(l3_p,l3_c,value(l3_c));
    l4_p0 = replace(l4_p,l4_c,value(l4_c));

    dx = [u_p1;u_p2;x1-x2; y1-y2];
    F2 = [
         sos(-dV_p0*dx - l3_p*reg)
         sos(l3_p)
         sos( dV_p0*dx + l4_p*reg)
         sos(l4_p)
         sum(abs(u_c1))>=1e-3 
         sum(abs(u_c2))>=1e-3
         sum(abs(u_c1))<=1 
         sum(abs(u_c2))<=1
         
         ];
    var2 = [u_c1;u_c2;l3_c;l4_c];

    [sol2,v,Q] = solvesos(F2,[],options,var2);
    sol2
    if sol2.problem == 0
        u_c10 = value(u_c1)'
        u_c20 = value(u_c2)'
        break;
    end
end

%%
delta = 0.5;
[x_mesh,y_mesh]= meshgrid(-5:delta:5,-5:delta:5);
ode_x = @(t,x) [
%                 x(1)*u_c10(1)+x(3)*u_c10(2)+x(2)*u_c10(3)+x(4)*u_c10(4)+...
%                 x(1)^2*u_c10(5)+x(1)*x(3)*u_c10(6)+x(3)^2*u_c10(7)+ ...
%                 x(1)*x(2)*u_c10(8)+x(3)*x(2)*u_c10(9)+x(2)^2*u_c10(10)+...
%                 x(1)*x(4)*u_c10(11)+x(3)*x(4)*u_c10(12)+x(2)*x(4)*u_c10(13)+...
%                 x(4)^2*u_c10(14);
%                 
%                 x(1)*u_c20(1)+x(3)*u_c20(2)+x(2)*u_c20(3)+x(4)*u_c20(4)+...
%                 x(1)^2*u_c20(5)+x(1)*x(3)*u_c20(6)+x(3)^2*u_c20(7)+ ...
%                 x(1)*x(2)*u_c20(8)+x(3)*x(2)*u_c20(9)+x(2)^2*u_c20(10)+...
%                 x(1)*x(4)*u_c20(11)+x(3)*x(4)*u_c20(12)+x(2)*x(4)*u_c20(13)+...
%                 x(4)^2*u_c20(14);
                
                u_c10*[x;1];u_c20*[x;1];

%                 u_c10*x;u_c20*x;
                x(1)-x(3);x(2)-x(4)];
[t,x_t] = ode45(ode_x,[0,1],10*rand(4,1));

figure;
plot(x_t(:,1),x_t(:,2),'r'); hold on;
plot(x_t(:,3),x_t(:,4),'b'); 
plot(x_t(1,1),x_t(1,2),'xr'); 
plot(x_t(1,3),x_t(1,4),'xb'); 
legend('evader','pursuer');
grid on;
axis equal;

x_p = [0;0];
for i = 1:length(x_mesh(:))
    temp = ode_x(0,[x_mesh(i);y_mesh(i);x_p(1);x_p(2)]);
    U_q(i,1) = temp(1);
    V_q(i,1) = temp(2);
end
figure;
quiver(x_mesh(:),y_mesh(:),U_q,V_q);
title(['Policy when the pursuer is at the origin'])
x_p'
% %%
% dX =[];
% for i = 1:length(x_mesh(:))
%     X(i,:) = [x_mesh(i), y_mesh(i)];
%     dX(i,:) = [u_c0(1)*x_mesh(i)+u_c0(2)*y_mesh(i); x_mesh(i)-y_mesh(i)];
% end 
% if length(V_c0) < 15
%     V_c0  = [0;V_c0];
% end
% V = V_c0(2)*X(:,1)+V_c0(3)*X(:,2)+V_c0(4)*X(:,1).^2+V_c0(5)*X(:,1).*X(:,2)+V_c0(6)*X(:,2).^2 ... ;
%     +X(:,1).^3*V_c0(7)+X(:,1).^2 .*X(:,2)*V_c0(8)+X(:,1).*X(:,2).^2*V_c0(9)+X(:,2).^3*V_c0(10)+ ...
%     X(:,1).^4*V_c0(11)+X(:,1).^3 .*X(:,2)*V_c0(12)+X(:,1).^2 .*X(:,2).^2*V_c0(13)+ ...
%     X(:,1).*X(:,2).^3*V_c0(14)+X(:,2).^4*V_c0(15);
% 
% dVdX(:,1) = V_c0(2) + 2*V_c0(4)*X(:,1) + V_c0(5)*X(:,2)  ...;
%        +3*X(:,1).^2*V_c0(7) + 2*X(:,1).*X(:,2)*V_c0(8) + X(:,2).^2*V_c0(9)+ ...
%        4*X(:,1).^3*V_c0(11) + 3*X(:,1).^2 .*X(:,2)*V_c0(12) + 2*X(:,1).*X(:,2).^2*V_c0(13)+ ...
%        X(:,2).^3*V_c0(14);
% dVdX(:,2) = V_c0(3) + V_c0(5)*X(:,1) + 2*V_c0(6)*X(:,2)  ...;
%        + X(:,1).^2 *V_c0(8)+ 2* X(:,1).*X(:,2)*V_c0(9)+3*X(:,2).^2*V_c0(10)+ ...
%        X(:,1).^3 *V_c0(12)+2*X(:,1).^2 .*X(:,2)*V_c0(13)+ ...
%        3*X(:,1).*X(:,2).^2*V_c0(14)+4*X(:,2).^3*V_c0(15);
% 
% 
% dV = sum(dVdX.*dX,2);
% 
% Map = x_mesh*0+1; Map(V<=0)=0;Map(dV>0) = 0;
% 
% V = reshape(V,size(x_mesh));
% dV = reshape(dV,size(x_mesh));
% 
% figure;
% surf(x_mesh,y_mesh,V);
% title('Lyapunov');
% axis([-5 5 -5 5])
% 
% figure;
% surf(x_mesh,y_mesh,dV);
% title('Grad Lyapunov');
% axis([-5 5 -5 5 ])
%  
% figure;
% surf(x_mesh,y_mesh,Map);
% title('Region with V>0 and dV<0');
% 
% if sol1.problem~=0
%     ep_system
% end
