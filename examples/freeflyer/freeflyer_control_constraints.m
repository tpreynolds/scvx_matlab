function f = freeflyer_control_constraints(u,auxdata)

% partial derivative
T = u(auxdata.id_T);
M = u(auxdata.id_M);

% norm(T) - T_max	<= 0.0;
% norm(M) - M_max	<= 0.0;

f = [ norm(T) - auxdata.T_nrm_max;
      norm(M) - auxdata.M_nrm_max ];

end

