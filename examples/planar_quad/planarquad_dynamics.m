function dx = planarquad_dynamics(input,t,x,u,p)

m = input.auxdata.m;
J = input.auxdata.J;
l = input.auxdata.l;

% r = x(1:2);
v = x(3:4);
a = x(5);
w = x(6);

% compute drag
[fD,~] = planarquad_drag(input,x);

F = [ -sin(a), -sin(a), -sin(a);
       cos(a),  cos(a),  cos(a) ];

% dynamics
dr = v;
dv = (1./m).*( F * u + fD ) + [ 0.0; input.auxdata.g];
da = w;
dw = (1/J).*dot([0;-l;l],u);

dx = [ dr; dv; da; dw ];

end
