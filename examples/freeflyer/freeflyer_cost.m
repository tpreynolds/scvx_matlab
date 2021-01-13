function J = freeflyer_cost(x,u,p,N)
scale = 1e0;
% J = 0.0; %norm( x(13*(N-1)+(1:3)) - [11.3; 6.0; 4.5;] );
J  = 0.0;
nu = numel(u)/N;
for k = 1:N-1
    uk  = u(nu*(k-1)+(1:nu));
    ukp = u(nu*k+(1:nu));
    J = J + 0.5 * scale * ( dot(uk,uk) + dot(ukp,ukp) );
end

end

