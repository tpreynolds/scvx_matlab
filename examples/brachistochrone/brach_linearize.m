function [A,B,C] = brach_linearize(input,t,x,u,p)

g = input.auxdata.g;

v = x(3);

A = zeros(input.nx,input.nx);
A(1,3) = sin(u);
A(2,3) = cos(u);

B = zeros(input.nx,input.nu);
B(1,1) =  v.*cos(u);
B(2,1) = -v.*sin(u);
B(3,1) = -g.*sin(u);

C = input.dynamics(input,t,x,u,p);

end

