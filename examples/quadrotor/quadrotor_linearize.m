function [A,B,C] = quadrotor_linearize(input,t,x,u,p)

id_r = input.auxdata.id_r;
id_v = input.auxdata.id_v;
id_w = input.auxdata.id_w;

A          = zeros(input.nx,input.nx);
A(id_r,id_v) = eye(3);

B          = zeros(input.nx,input.nu);
B(id_v,id_w) = eye(3);

C = input.dynamics(input,t,x,u,p);

end

