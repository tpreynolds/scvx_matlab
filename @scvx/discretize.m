function discretize(obj)
%DISCRETIZE
%
% T. Reynolds -- RAIN Lab

    nx  = obj.nx;
    nu  = obj.nu;
    np  = obj.np;
    N   = obj.ctrl.N;
    tau = obj.auxdata.tau;

    % initialize data
    EH = scvx.set_block(zeros(nx*N,nx*N),eye(nx),1,1);
    BE = zeros(nx*N,nu*N);
    ES = zeros(nx*N,np);
    AR = zeros(nx*N,1);
    feas        = true;     % default value
    first_iter  = false;    % default value
    if (isempty(obj.output))
        first_iter = true;
    end
    
    % compute initial condition
    if (first_iter)
        x0 = obj.iguess.x(:,1);
        p  = obj.iguess.p;
    else
        x0 = scvx.get_block(obj.output.x,1,1,N,1);
        p  = obj.output.p;
    end
    A0 = eye(nx);
    B0 = zeros(2*nx*nu,1);
    S0 = zeros(nx,np);
    R0 = zeros(nx,1);
    P0 = [ x0(:); A0(:); B0(:); S0(:); R0(:) ];
    
    for k = 1:N-1
       
        tspan = linspace(tau(k),tau(k+1),obj.ctrl.Nsub);
        if (first_iter)
            u = obj.iguess.u(:,k:k+1);
        else
            u = [ scvx.get_block(obj.output.u,k,1,N,1), ...
                  scvx.get_block(obj.output.u,k+1,1,N,1) ];
        end
        % integrate 
        P = scvx.rk4(@(t,Z)deriv(t,Z,u,p),tspan,P0);
        % grab the current state at node k+1
        if (first_iter)
            xkp = obj.iguess.x(:,k+1);
        else
            xkp = scvx.get_block(obj.output.x,k+1,1,N,1);
        end
        % extract elements of stacked vector P
        xP  = P(end,1:nx);
        AP  = P(end,nx+(1:nx*nx));
        BmP = P(end,nx*(nx+1)+(1:nx*nu));
        BpP = P(end,nx*(nx+nu+1)+(1:nx*nu));
        SP  = P(end,nx*(nx+2*nu+1)+(1:nx*np));
        RP  = P(end,nx*(nx+2*nu+np+1)+(1:nx));
        % reshape to matrices
        xd  = reshape(xP,nx,1);
        Ad  = reshape(AP,nx,nx);
        Bmd = Ad * reshape(BmP,nx,nu);
        Bpd = Ad * reshape(BpP,nx,nu);
        Sd  = Ad * reshape(SP,nx,np);
        Rd  = Ad * reshape(RP,nx,1);
        % fill up concatenated matrices
        EH = scvx.set_block(EH,Ad,k+1,k);
        BE = scvx.set_block(BE,Bmd,k+1,k);
        BE = scvx.set_block(BE,Bpd,k+1,k+1);
        ES = scvx.set_block(ES,Sd,k+1,1);
        AR = scvx.set_block(AR,Rd,k+1,1);
        % check residual and assess feasibility
        defect = norm(xkp-xd,2);
        if (defect>obj.ctrl.feas_tol && feas)
            feas = false;
        end
        obj.defects(k) = defect;
        % reset initial condition
        P0(1:nx) = xkp;
    end

    obj.output.Ad   = EH;
    obj.output.Bd   = BE;
    obj.output.Sd   = ES;
    obj.output.Rd   = AR;
    obj.output.feas = feas;
    
    function dZ = deriv(tau,Z,u,p)
        
        x_tau = Z(1:nx);
        u_tau = interp1([tspan(1), tspan(end)],u',tau,'linear')';
        PHI   = reshape(Z(nx+(1:nx*nx)),nx,nx);
        lm    = (tspan(end)-tau)/(tspan(end)-tspan(1));
        lp    = (tau-tspan(1))/(tspan(end)-tspan(1));
        
        f = obj.dynamics(obj,tau,x_tau,u_tau,p);
        [A,B,C] = obj.linearize(obj,tau,x_tau,u_tau,p);
        
        Bm   = lm * B;
        Bp   = lp * B;
        R    = -A*x_tau-B*u_tau;
        iPHI = eye(size(PHI))/PHI;
        
        if (obj.final_time_free)
            tf = p(1);
        else
            tf = 1.0;
        end
        dx = tf * f;
        dPHI = (tf * A) * PHI;
        dBm  = iPHI * (tf * Bm);
        dBp  = iPHI * (tf * Bp);
        dS   = iPHI * C;
        dR   = iPHI * (tf * R);
        
        dZ = [ dx; dPHI(:); dBm(:); dBp(:); dS(:); dR(:) ];
        
    end

end
