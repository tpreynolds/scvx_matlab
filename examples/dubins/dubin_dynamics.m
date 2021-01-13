function dx = dubin_dynamics(input,t,x,u,p)

a = x(3);
v = u(1);
w = u(2);

dxx = v.*sin(a);
dyy = v.*cos(a);
daa = w;

dx = [ dxx; dyy; daa ];

end

