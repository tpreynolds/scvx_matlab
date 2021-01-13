function dx = planarlanding_dynamics(input,t,x,u,p)
%PLANARLANDING_DYNAMICS     
%
% Continuous dynamic equations for the planar landing problem.
%
% Syntax is dx = obj.dynamics(obj,t,x,u,p) where obj is an scvx object.
%
% T. Reynolds -- RAIN Lab

m = x(1);
r = x(2:3);
v = x(4:5);
a = x(6);
w = x(7);

G = u(1);
T = u(2);

d = [ -sin(a); cos(a) ];

% dynamics
dm = -input.auxdata.alpha .* G;
dr = v;
dv = (G./m).*d + [ 0.0; input.auxdata.g];
da = w;
dw = T./input.auxdata.inertia;

dx = [ dm; dr; dv; da; dw ];

end
