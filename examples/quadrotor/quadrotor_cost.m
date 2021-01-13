function J = quadrotor_cost(x,u,p,N)
scale = 1e2;
J  = 0.0;
nu = numel(u)/N;
for k = 1:N-1
    uk  = u(nu*(k-1)+(1:nu));
    ukp = u(nu*k+(1:nu));
    J = J + 0.5 * scale * ( norm(uk(4),2) + norm(ukp(4),2) );
end

end

