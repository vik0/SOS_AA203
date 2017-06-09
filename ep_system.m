close all
clear 
clc

sdpvar x1 x2
x = [x1; x2];

[V_p,V_c] = polynomial(x,4,1);
dV_p = jacobian(V_p,x);

r0 = 0.5;
reg = (x1-x2)^2 - r0^2;

[u_p,u_c] = polynomial(x,1,1);

[l1_p,l1_c] = polynomial(x,4);
[l2_p,l2_c] = polynomial(x,4);
[l3_p,l3_c] = polynomial(x,4);
[l4_p,l4_c] = polynomial(x,4);
u_c0 = randn(size(u_c));

u_p0 = replace(u_p,u_c,u_c0);
dx0 = [u_p0; x1-x2];
F1 = [
    sos(V_p - l1_p*reg)
    sos(l1_p)
    sos(-V_p + l2_p*reg)
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
sdpvar a b
dx = [a*x1 + b*x2;x1-x2];
F2 = [
     sos(-dV_p0*dx - l3_p0*reg)
     sos( dV_p0*dx + l4_p0*reg)
     ];
var2 = [a;b];

[sol2,v,Q] = solvesos(F2,[],options,var2);
sol2
if sol2.problem == 0
    u_c0 = value([a;b]')
end

%%
delta = 0.5;
[x_mesh,y_mesh]= meshgrid(-5:delta:5,-5:delta:5);
figure;
uv(:,1) = u_c0(1)*x_mesh(:)+u_c0(2)*y_mesh(:);
uv(:,2) = x_mesh(:)-y_mesh(:);
quiver(x_mesh(:),y_mesh(:),uv(:,1),uv(:,2)); hold on;
plot([-5,5],[-5,5],'r');
title('quiver plot');
axis equal;


dX =[];
for i = 1:length(x_mesh(:))
    X(i,:) = [x_mesh(i), y_mesh(i)];
    dX(i,:) = [u_c0(1)*x_mesh(i)+u_c0(2)*y_mesh(i); x_mesh(i)-y_mesh(i)];
end 
if length(V_c0) < 15
    V_c0  = [0;V_c0];
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

Map = x_mesh*0+1; Map(V<=0)=0;Map(dV>0) = 0;

V = reshape(V,size(x_mesh));
dV = reshape(dV,size(x_mesh));

figure;
surf(x_mesh,y_mesh,V);
title('Lyapunov');
axis([-5 5 -5 5])

figure;
surf(x_mesh,y_mesh,dV);
title('Grad Lyapunov');
axis([-5 5 -5 5 ])
 
figure;
surf(x_mesh,y_mesh,Map);
title('Region with V>0 and dV<0');

if sol1.problem~=0
    ep_system
end
