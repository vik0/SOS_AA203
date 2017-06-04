close all
clear 
clc

sdpvar x1 y1 x2 y2
x = [x1;y1;x2;y2];

v1max = 1; v2max = 0.1;

A = [     -1,  0, 0, 0;
           0, -1, 0, 0;
           1,  0,-1, 0;
           0,  1, 0,-1];

dx_1 = A*x;
dx(1:2,1) = v1max* dx_1(1:2,1)/norm(dx_1(1:2,1));
dx(3:4,1) = v2max* dx_1(3:4,1)/norm(dx_1(3:4,1));

[V_p, V_c]  = polynomial(x,2,1);
[l1_p,l1_c] = polynomial(x,2);
[l2_p,l2_c] = polynomial(x,2);
dV_p = jacobian(V_p,x);

% safe region 
R0 = 0.1;
safe = R0- x2^2 -y2^2;


options = sdpsettings('sos.newton',1,'sos.congruence',1,'verbose',0);

F = [sos(V_p + l1_p*safe)
     sos(l1_p)
     sos(-dV_p*dx + l2_p*safe)
     sos(l2_p)
     sum(V_c)~=0];
 
if (value(x1^2 + y1^2)>R0)
    disp('fi')
    F = [F;
     % constraint that the distance should be greater than Ro
     ((x1-x2)^2 + (y1-y2)^2 -R0 )>=0
     ((x2^2+y2^2)-a*(x1^2+y1^2))>0
     ];
end

var = [V_c; l1_c; l2_c];

[sol,v,Q] = solvesos(F,[],[],var);
sol
% var(1)+x1*var(2)+y1*var(3)+x2*var(4)+y2*var(5)+x1^2*var(6)+x1*y1*var(7)+y1^2*var(8)+x1*x2*var(9)+y1*x2*var(10)+x2^2*var(11)+x1*y2*var(12)+y1*y2*var(13)+x2*y2*var(14)+y2^2*var(15)
value(V_c)'
