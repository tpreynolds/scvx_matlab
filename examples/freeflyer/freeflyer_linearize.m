function [A,B,C] = freeflyer_linearize(input,t,x,u,p)

id_r = input.auxdata.id_r;
id_v = input.auxdata.id_v;
id_q = input.auxdata.id_q;
id_w = input.auxdata.id_w;
id_T = input.auxdata.id_T;
id_M = input.auxdata.id_M;

q = x(id_q);
w = x(id_w);

m = input.auxdata.mass;
J = input.auxdata.inertia;

dfq_dq = 0.5 * [  0.0,  w(3), -w(2),  w(1);
                -w(3),   0.0,  w(1),  w(2);
                 w(2), -w(1),   0.0,  w(3);
                -w(1), -w(2), -w(3),  0.0 ];
dfq_dw = 0.5 * [ q(4), -q(3),  q(2);
                 q(3),  q(4), -q(1);
                -q(2),  q(1),  q(4);
                -q(1), -q(2), -q(3) ];
dfw_dw = -J\( skew(w)*J - skew(J*w) );


A = zeros(input.nx,input.nx);
A(id_r,id_v) = eye(3);
A(id_q,id_q) = dfq_dq;
A(id_q,id_w) = dfq_dw;
A(id_w,id_w) = dfw_dw;

B = zeros(input.nx,input.nu);
B(id_v,id_T) = (1/m)*eye(3);
B(id_w,id_M) = J\eye(3);

C = input.dynamics(input,t,x,u,p);

end

