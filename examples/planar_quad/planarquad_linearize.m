function [A,B,C] = planarquad_linearize(input,t,x,u,p)

m = input.auxdata.m;
J = input.auxdata.J;
l = input.auxdata.l;

[~,dfD_dx] = planarquad_drag(input,x);

sa = sin(x(5));
ca = cos(x(5));

F = [ -sa, -sa, -sa;
       ca,  ca,  ca ];
dFda = [ -ca, -ca, -ca;
         -sa, -sa, -sa ];

A           = zeros(input.nx,input.nx);
A(1:2,3:4)  = eye(2);
A(3:4,:)    = dfD_dx;
A(3:4,5)    = (1/m) .* dFda * u;
A(5,6)      = 1.0;

B        = zeros(input.nx,input.nu);
B(3:4,:) = (1./m) .* F;
B(6,:)   = (1/J) .* [ 0, -l, l ];

C = input.dynamics(input,t,x,u,p);

end

