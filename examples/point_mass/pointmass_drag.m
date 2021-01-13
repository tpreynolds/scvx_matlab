function [fD,dfD_dx] = pointmass_drag(input,x)

q = - 0.5 * input.auxdata.rho * input.auxdata.Sd * input.auxdata.cd;

m = input.auxdata.m;
v = x(4:6);
speed = norm(v);

fD = q * speed .* v;

dfD_dx = zeros(3,input.nx);
if (speed>1e-12)
    dfD_dx(:,4:6) = (q/m) .* ( speed * eye(3) + (v*v')/speed );
end

end

