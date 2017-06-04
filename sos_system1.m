close all
clear all
clc

sdpvar x1 y1 x2 y2 t
[p1,c1] = polynomial([x1,y1],2,1);
[p2,c2] = polynomial([x2,y2],2,1);
j1 = jacobian(p1,[x1,y1]);
j2 = jacobian(p2,[x2,y2]);

A = -eye(2);
dxy1 = A*[x1;y1];
% dxy1 = dxy1/(norm(dxy1)+0.01);
dxy2 = [x1-x2; y1-y2];
% dxy2 = dxy2/(norm(dxy2)+0.01);

R0 = 1e-3;

F = [sos(p1)
     sos(-j1*dxy1)
     sos(p2)
     sos(-j2*dxy2)];
 
% if (value(x1^2 + y1^2)>R0)
%     F = [F;
%      % constraint that the distance should be greater than Ro
%      ((x1-x2)^2 + (y1-y2)^2 -R0 )>=0
%      ((x2^2+y2^2)-(x1^2+y1^2))>0
%      ];
% end

sol = solvesos(F,[],[],[c1;c2]);

sol
value(c1)'
value(c2)'