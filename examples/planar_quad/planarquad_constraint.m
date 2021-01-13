function [f,A] = planarquad_constraint(x)

% define size/shape
H  = diag([1,1]);
% define center
xc = [1;8];

diff = x(1:2) - xc;

f = 1 - norm(H * diff); % <= 0

% partial derivative
A = zeros(1,numel(x));
A(1:2) = - (diff'*(H'*H))/(norm(H*diff));

end

