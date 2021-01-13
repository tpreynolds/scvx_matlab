function dx = pointmass_dynamics(input,t,x,u,p)

m = input.auxdata.m;
v = x(4:6);

% compute drag
[fD,~] = pointmass_drag(input,x);

% dynamics
dr = v;
dv = (1./m).*(u + fD) + [ 0.0; 0.0; input.auxdata.g];

dx = [ dr; dv ];

end
