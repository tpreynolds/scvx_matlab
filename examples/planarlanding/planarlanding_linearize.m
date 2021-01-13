function [A,B,C] = planarlanding_linearize(input,t,x,u,p)
%LINEARIZE    
%
% Returns linearized matrices at point (t,x,u,p) for the planar landing
% problem.
%
% Syntax is [A,B,C] = obj.linearize(obj,t,x,u,p) where obj is an scvx
% object and outputs are
%   - A : Jacobian dfdx
%   - B : Jacobian dfdu
%   - C : Jacobian dfdp
%
% T. Reynolds -- RAIN Lab

m = x(1);
a = x(6);
G = u(1);

d   = [ -sin(a); cos(a) ];
da  = [ -cos(a); -sin(a) ];

A           = zeros(input.nx,input.nx);
A(4:5,1)    = -(G./m^2).*d; 
A(2:3,4:5)  = eye(2);
A(4:5,6)    = (G./m).*da;
A(6,7)      = 1.0;

B        = zeros(input.nx,input.nu);
B(1,1)   = -input.auxdata.alpha;
B(4:5,1) = (1./m).*d;
B(7,2)   = 1./input.auxdata.inertia;

C = input.dynamics(input,t,x,u,p);

end

