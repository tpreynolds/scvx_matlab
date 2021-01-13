function [A,B,C] = pointmass_linearize(input,t,x,u,p)

m = input.auxdata.m;

[~,dfD_dx] = pointmass_drag(input,x);

A           = zeros(input.nx,input.nx);
A(1:3,4:6)  = eye(3);
A(4:6,:)    = dfD_dx;

B        = zeros(input.nx,input.nu);
B(4:6,:) = (1./m)*eye(3);

C = input.dynamics(input,t,x,u,p);

end

