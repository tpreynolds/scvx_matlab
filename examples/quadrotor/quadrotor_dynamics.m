function dx = quadrotor_dynamics(input,t,x,u,p)

v = x(input.auxdata.id_v);
w = u(input.auxdata.id_w);

% dynamics
dr = v;
dv = w + [ 0.0; 0.0; input.auxdata.g];

dx = [ dr; dv ];

end
