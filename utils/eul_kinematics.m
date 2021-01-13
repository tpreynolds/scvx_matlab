function [T,dT_deul] = eul_kinematics(eul)

roll    = eul(1);
pitch   = eul(2);

sr = sin(roll);
cr = cos(roll);
tp = tan(pitch);
sep = sec(pitch);

T = [ 1.0,  sr*tp,   cr*tp;
      0.0,  cr,     -sr;
      0.0,  sr*sep,  cr*sep ];
  
dT_deul = zeros(3,3,3);

dT_deul(:,:,2) = [ cr*tp,   sr*sep*sep, 0.0;
                   -sr,     0.0,        0.0;
                   cr*sep,  sr*tp*sep,  0.0 ];
               
dT_deul(:,:,3) = [ -sr*tp,  cr*sep*sep, 0.0;
                   -cr,     0.0,        0.0;
                   -sr*sep, cr*tp*sep,  0.0 ];

end

