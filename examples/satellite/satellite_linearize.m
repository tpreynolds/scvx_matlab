function [A,B,C] = satellite_linearize(input,t,x,u,p)

J = input.auxdata.J;
id_a = input.auxdata.id_a;
id_w = input.auxdata.id_w;

eul = x(id_a);
wB  = x(id_w);

[T,dT_deul] = eul_kinematics(eul);

A = zeros(input.nx,input.nx);
A(id_a,id_a) = dT_deul(:,:,1) * wB(1) ...
                       + dT_deul(:,:,2) * wB(2) ...
                       + dT_deul(:,:,3) * wB(3);
A(id_a,id_w) = T;
A(id_w,id_w) = J\( skew(J*wB) - skew(wB)*J );

B = zeros(input.nx,input.nu);
B(id_w,:) = J\eye(3);

C = input.dynamics(input,t,x,u,p);

end

