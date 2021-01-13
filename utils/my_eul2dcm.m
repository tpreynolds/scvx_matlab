function [C_I2B,dC_I2B] = my_eul2dcm(Theta)

cr = cos(Theta(1));
sr = sin(Theta(1));
cp = cos(Theta(2));
sp = sin(Theta(2));
cy = cos(Theta(3));
sy = sin(Theta(3));

C_I2B = [ cp*cy,          cp*sy,          -sp;
          sr*sp*cy-cr*sy, sr*sp*sy+cr*cy, sr*cp;
          cr*sp*cy+sr*sy, cr*sp*sy-sr*cy, cr*cp ];
      
if (nargout>1)
   dC_I2B(:,:,1) = [ 0.0, 0.0, 0.0;
                     cr*sp*cy+sr*sy, cr*sp*sy-sr*cy, cr*cp;
                     -sr*sp*cy+cr*sy, -sr*sp*sy-cr*cy, -sr*cp ];
   dC_I2B(:,:,2) = [ -sp*cy, -sp*sy, -cp;
                    sr*cp*cy, sr*cp*sy, -sr*sp;
                   -cr*cp*cy, cr*cp*sy, -cr*sp ];
   dC_I2B(:,:,3) = [ -cp*sy, cp*cy, 0.0;
                     -sr*sp*sy-cr*cy, sr*sp*cy-cr*sy, 0.0;
                     -cr*sp*sy+sr*cy, cr*sp*cy+sr*sy, 0.0 ];
end

end

