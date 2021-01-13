function varxu = solve_socp(obj,scaling,iter)
%SOLVE_SOCP
%
% T. Reynolds

nx = obj.nx;
nu = obj.nu;
np = obj.np;
N  = obj.ctrl.N;
first   = 1:nx;
last    = nx*(N-1)+1:nx*N;
isptr   = strcmp(obj.ctrl.algo,'ptr');

% reference values
wvc     = obj.ctrl.wvc;
wtr     = obj.ctrl.wtr;
wtrp    = obj.ctrl.wtrp;
if (~isptr)
    tr = obj.output.tr;
end
if (iter < 2) 
    xref = obj.iguess.x(:);
    uref = obj.iguess.u(:);
    pref = obj.iguess.p;
else
    xref = obj.output.x;
    uref = obj.output.u;
    pref = obj.output.p;
end
EH = obj.output.Ad;
BE = obj.output.Bd;
ES = obj.output.Sd;
AR = obj.output.Rd;

% scaling matrices
Sx = scaling.Sx;
cx = scaling.cx;
SX = scaling.SX;
cX = scaling.cX;
Su = scaling.Su;
cu = scaling.cu;
SU = scaling.SU;
cU = scaling.cU;
Sp = scaling.Sp;
cp = scaling.cp;
iSx = scaling.iSx;
iSu = scaling.iSu;
iSp = scaling.iSp;

% constraint things
x_min = obj.bnds.x_min;
x_max = obj.bnds.x_max;
u_min = obj.bnds.u_min;
u_max = obj.bnds.u_max;
nx_constraints_path = numel(obj.bnds.path.x);
nu_constraints_path = numel(obj.bnds.path.u);

cvx_tic;
cvx_begin quiet
    cvx_solver(obj.ctrl.solver)
    cvx_precision('low')
    
    variable xb(nx*N,1)
    variable ub(nu*N,1)
    variable pb(np,1)   % even in np=0, this is just an empty vector
    variable vc(nx*N,1)
    
    cost = obj.cost(xb,ub,pb,N) + wvc * norm(vc,1);
    if (np>0)
        if (isptr)
            variable tr(N) nonnegative
            variable tr_p nonnegative
            cost = cost + wtr*norm(tr,1) + wtrp*norm(tr_p,1);
        end
    else % no parameters
        if (isptr)
            variable tr(N) nonnegative
            cost = cost + wtr*norm(tr,1);
        end
    end
    minimize( cost )
    
    subject to
    
    % initial conditions
    (Sx*xb(first)+cx) <= obj.bnds.init.x_max;
    (Sx*xb(first)+cx) >= obj.bnds.init.x_min;
    
    % final conditions
    (Sx*xb(last)+cx) <= obj.bnds.trgt.x_max;
    (Sx*xb(last)+cx) >= obj.bnds.trgt.x_min;
    
    % bounds on parameters
    (Sp*pb+cp) <= obj.bnds.p_max;
    (Sp*pb+cp) >= obj.bnds.p_min;
    
    % dynamics
    (SX*xb+cX) == EH*(SX*xb+cX) + BE*(SU*ub+cU) + ES*(Sp*pb+cp) + AR + vc;
    
    % parameter trust region
    if (np>0 && isptr)
        (pb-(iSp*(pref-cp)))'*(pb-(iSp*(pref-cp))) <= tr_p;
    end
    
    % time loop
    for k = 1:N
       
        xbk   = xb(nx*(k-1)+1:nx*k);
        xrefk = xref(nx*(k-1)+1:nx*k);
        ubk   = ub(nu*(k-1)+1:nu*k);
        urefk = uref(nu*(k-1)+1:nu*k);
        
        % trust region
        if (isptr)
            % quadratic trust region
            (xbk - (iSx*(xrefk-cx)))'*(xbk-(iSx*(xrefk-cx))) ...
                + (ubk-(iSu*(urefk-cu)))'*(ubk-(iSu*(urefk-cu))) <= tr(k);
            % linear trust region
%             -tr(k) <= xbk - (iSx*(xrefk-cx)) <= tr(k);
%             -tr(k) <= ubk - (iSu*(urefk-cu)) <= tr(k);
        else
            % quadratic trust region
%             (xbk - (iSx*(xrefk-cx)))'*(xbk-(iSx*(xrefk-cx))) ...
%                + (ubk-(iSu*(urefk-cu)))'*(ubk-(iSu*(urefk-cu))) <= tr;
           % quadratic trust region (state only)
%             (xbk - (iSx*(xrefk-cx)))'*(xbk-(iSx*(xrefk-cx))) <= tr;
            % linear trust region
            -tr <= xbk - (iSx*(xrefk-cx)) <= tr;
            -tr <= ubk - (iSu*(urefk-cu)) <= tr;
        end
        
        % linear bounds
        (Sx*xbk+cx) <= x_max;
        (Sx*xbk+cx) >= x_min;
        (Su*ubk+cu) <= u_max;
        (Su*ubk+cu) >= u_min;
        
        % state path constraints
        for xcnstr = 1:nx_constraints_path
            if (obj.bnds.path.x_cvx{xcnstr})
                obj.bnds.path.x{xcnstr}( Sx*xbk+cx ) <= 0.0;
            else
                TS = obj.bnds.path.x_lin{xcnstr};
                TS{1}(k) + TS{2}(k,:) * (Sx*xbk+cx) <= 0.0;
            end
        end
                       
        % control path constraints
        for ucnstr = 1:nu_constraints_path
            if (obj.bnds.path.u_cvx{ucnstr})
                obj.bnds.path.u{ucnstr}( Su*ubk+cu ) <= 0.0;
            else
                TS = obj.output.bnds.u{ucnstr};
                TS{1} + TS{2} * (Su*ubk+cu) <= 0.0;
            end
        end
        
    end
cvx_end
timing = cvx_toc;
    
obj.output.x    = SX*xb+cX;
obj.output.u    = SU*ub+cU;
obj.output.p    = Sp*pb+cp;
obj.output.vc   = vc;
obj.output.tr   = tr;
obj.output.cost = cvx_optval;

% compute max (scaled) change in state and/or control
temp = 0.0;
for k = 1:N
   tempx = norm(xb(nx*(k-1)+(1:nx)) - iSx*(xref(nx*(k-1)+(1:nx))-cx),inf);
   tempu = norm(ub(nu*(k-1)+(1:nu)) - iSu*(uref(nu*(k-1)+(1:nu))-cu),inf);
   temp  = max([temp,tempx,tempu]);    
end
varxu = temp;

% save some data from the solution
if (iter>1)
    obj.output.data.cum_time = obj.output.data.cum_time + sum(timing(4:5));
else
    obj.output.data.cum_time = sum(timing(4:5));
end
switch cvx_status
    case 'Solved'
        obj.output.data.status{iter} = 0;
    case 'Inaccurate/Solved'
        obj.output.data.status{iter} = 1;
    case {'Failed','Unbounded','Inaccurate/Unbounded'}
        obj.output.data.status{iter} = 2;
    case {'Infeasible','Inaccurate/Infeasible'}
        obj.output.data.status{iter} = 3;
    otherwise
        obj.output.data.status{iter} = 4;
end
obj.output.data.slvitr = cvx_slvitr;

end

