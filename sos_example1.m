% example for checking if polynomial can be written as SOS

close all
clear all
clc

x = sdpvar(1,1);
y = sdpvar(1,1);
z = sdpvar(1,1);

p = 12+y^2-2*x^3*y+2*y*z^2+x^6-2*x^3*z^2+z^4+x^2*y^2; 
% p = x^2 + y^2;

options = sdpsettings('sos.newton',1,'sos.congruence',1,'sos.numblkdg',1e-4);
[sol,v,Q] = solvesos(sos(p),[],options);

disp('PSD Matrix for SOS formulation is ');
Q{1}
disp('Monomials used are ');
sdisplay(v{1}')

disp('succesful execution suggests that polynomail can written as SOS');
disp('Experiment by tweaking coefficients of the polynomial');

