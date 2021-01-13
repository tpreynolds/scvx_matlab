function [A,B,C] = dubin_linearize(input,t,x,u,p)

a = x(3);
v = u(1);

A = zeros(input.nx,input.nx);
A(1,3) = v*cos(a);
A(2,3) = -v*sin(a);

B = zeros(input.nx,input.nu);
B(1,1) = sin(a);
B(2,1) = cos(a);
B(3,2) = 1.0;

C = input.dynamics(input,t,x,u,p);

end

