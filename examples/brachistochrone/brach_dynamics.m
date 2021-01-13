function dx = brach_dynamics(input,t,x,u,p)

g = input.auxdata.g;

v = x(3);

dx = v.*sin(u);
dy = v.*cos(u);
dv = g.*cos(u);

dx = [ dx; dy; dv ];

end

