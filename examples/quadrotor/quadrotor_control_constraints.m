function f = quadrotor_control_constraints(u,auxdata)

% partial derivative
w = u(1:3);
G = u(4);

% norm(w) - G                     <= 0.0;
% G - auxdata.u_nrm_max           <= 0.0;
% auxdata.u_nrm_min - G           <= 0.0; 
% G*cos(auxdata.tilt_max) - w(3)  <= 0.0;

f = [ norm(w) - G;
      G - auxdata.u_nrm_max;
      auxdata.u_nrm_min - G;
      G*cos(auxdata.tilt_max)-w(3) ];

end

