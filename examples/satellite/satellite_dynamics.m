function dx = satellite_dynamics(input,t,x,u,p)

J = input.auxdata.J;

eul = x(input.auxdata.id_a);
wB  = x(input.auxdata.id_w);

T = eul_kinematics(eul);

deul = T * wB;
dwB  = J\( u - skew(wB)*J*wB );

dx = [ deul; dwB ];

end

