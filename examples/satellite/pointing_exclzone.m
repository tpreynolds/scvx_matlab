function [f,dfdx] = pointing_exclzone(x)

theta = deg2rad(40);
xI = [1;0;0];
yB = [0;1;0];

% A = xI*yB' + yB*xI' - (dot(xI,yB) + cos(theta)) * eye(3);
% b = cross(xI,yB);
% c = dot(xI,yB) - cos(theta);
% M = [ A,b;b',c ] + 2*eye(4);

[C_I2B,dC_I2B] = my_eul2dcm(x(1:3));

f = xI' * C_I2B' * yB - cos(theta);
if (nargout>1)
    dfdx = zeros(1,numel(x));
    dfdx(1:3) = xI' * ( dC_I2B(:,:,1)'*yB(1) ...
                      + dC_I2B(:,:,2)'*yB(2) ...
                      + dC_I2B(:,:,3)'*yB(3) ); 
      
end

