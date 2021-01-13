function C = rot(a,type)
switch type
    case 'x'
        C = [ 1 0 0;
              0 cosd(a) sind(a);
              0 -sind(a) cosd(a) ];
    case 'y'
        C = [ cosd(a) 0 -sind(a);
                0   1   0;
                sind(a)  0   cosd(a) ];
    case 'z'
        C = [ cosd(a) sind(a) 0;
              -sind(a) cosd(a) 0;
              0     0    1 ];
    otherwise
        error('ROT: input type must be one of x, y, or z')
end
end
