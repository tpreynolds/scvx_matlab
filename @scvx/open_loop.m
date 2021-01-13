function result = open_loop(obj)
%OPEN_LOOP

nx = obj.nx;
N  = obj.ctrl.N;
result = struct;

% number of points in the integration 
Nfull = 3 * obj.ctrl.Nsub * N;

% initial condition
x0 = obj.output.x(1:nx);

% time span
if (obj.initial_time_free && obj.final_time_free)
    t0 = obj.output.p(2);
    tf = obj.output.p(1);
elseif (obj.final_time_free)
    t0 = obj.bnds.init.t_min;
    tf = obj.output.p(1);
elseif (obj.initial_time_free)
    t0 = obj.output.p(1);
    tf = obj.bnds.trgt.t_min;
end
tspan = linspace(t0,tf,Nfull);

% get control matrix (reshape to nu by N) and parameters p
u   = reshape( obj.output.u, obj.nu, N );
ut  = linspace(t0,tf,N);
p   = obj.output.p;

% propagate dynamics
Z = scvx.rk4(@(t,x)deriv(t,x,u,p),tspan,x0);

% compute final state error
xf  = obj.output.x(nx*(N-1)+(1:nx));
Zf  = Z(end,:)';
err = norm(xf-Zf,2);

% output of the open loop propagation
result.x        = Z';
result.t        = tspan;
result.err      = err;
result.Nfull    = Nfull;

    function dx = deriv(t,x,u,p)
        
        % interpolate control
        u_t = interp1(ut,u',t,'linear')';
        
        % query dynamics
        dx = obj.dynamics(obj,t,x,u_t,p);
        
    end

end

