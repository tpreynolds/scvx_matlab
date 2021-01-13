function J = planarquad_cost(x,u,p,N)
scale = 1e1;
J  = 0.0;
nu = numel(u)/N;
for k = 1:N-1
    uk  = u(nu*(k-1)+(1:nu));
    ukp = u(nu*k+(1:nu));
    J = 0.5 * scale * ( norm(uk,1) + norm(ukp,1) );
end

end

