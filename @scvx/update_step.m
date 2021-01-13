function [reject,change] = update_step(obj)
%UPDATE_STEP
%
% [reject,change] = update_step(obj) performs the trust region update step
% for the SCVX algorithm. 
% 
%   The boolean reject indicates if the step was accepted or rejected. 
% If accepted, obj.output is written into obj.best_output. If rejected, 
% obj.best_output is read back into obj.output and the trust region is 
% shrunk. 
% 
% The output change is a char that indicates what happened with the trust
% region. Possible values are:
%   S : trust region was shrunk
%   K : trust region was kept
%   G : trust region was grown
%   L : trust region was shrunk but is at its lower bound
%   U : trust region was grown but is at its upper bound
%
% T. Reynolds -- RAIN Lab

% compute new costs
cost  = obj.cost(obj.output.x,obj.output.u,obj.output.p,obj.ctrl.N);
new_L = cost + obj.ctrl.wvc * norm(obj.output.vc,1);
new_J = cost + obj.ctrl.wvc * sum(obj.defects);

dL = obj.last_J - new_L;
dJ = obj.last_J - new_J;

% compute performance metric
rho = dJ/dL;

% update trust region radius
if ( rho < obj.ctrl.rho0 )
    reject = true;
    obj.output.tr = obj.output.tr / obj.ctrl.alpha;
    change = 'S';
else
    reject = false;
    obj.last_J = new_J;
    if ( rho < obj.ctrl.rho1 )
        % shrink
        obj.output.tr = obj.output.tr / obj.ctrl.alpha;
        change = 'S';
    elseif ( (rho >= obj.ctrl.rho1) && (rho < obj.ctrl.rho2) )
        change = 'K';
    else
        % grow
        obj.output.tr = obj.output.tr * obj.ctrl.beta;
        change = 'G';
    end
end

% saturate trust region to lower/upper bound 
if (obj.output.tr < obj.ctrl.tr_lb)
    obj.output.tr = obj.ctrl.tr_lb;
    change = 'L';
elseif (obj.output.tr > obj.ctrl.tr_ub)
    obj.output.tr = obj.ctrl.tr_ub;
    change = 'U';
end

% if rejected, recover the last best output. If accepted, set the new best
% output
if (reject)
    obj.get_best_output();
else
    obj.set_best_output();
end
    
end

