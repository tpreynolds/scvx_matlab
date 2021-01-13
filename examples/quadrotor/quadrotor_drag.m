function [fD,dfD_dx] = planarquad_drag(input,x)

q = - 0.5 * input.auxdata.rho * input.auxdata.Sd * input.auxdata.cd;

m = input.auxdata.m;
v = x(3:4);
speed = norm(v);

fD = q * speed .* v;

dfD_dx = zeros(2,input.nx);
if (speed>1e-12)
    dfD_dx(:,3:4) = (q/m) .* ( speed * eye(2) + (v*v')/speed );
end

end

