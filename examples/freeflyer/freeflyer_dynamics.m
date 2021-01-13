function dx = freeflyer_dynamics(input,t,x,u,p)

v = x(input.auxdata.id_v);
q = x(input.auxdata.id_q);
w = x(input.auxdata.id_w);
T = u(input.auxdata.id_T);
M = u(input.auxdata.id_M);

m = input.auxdata.mass;
J = input.auxdata.inertia;

% dynamics
dr = v;
dv = T/m;
dq = 0.5 * [ q(4)*eye(3)+skew(q(1:3)), q(1:3); -q(1:3)', q(4) ] * [w;0];
dw = J\( M - skew(w) * J * w);

dx = [ dr; dv; dq; dw ];

end
