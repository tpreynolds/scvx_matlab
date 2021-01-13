function varargout = solve(obj)
%SOLVE
%
% Main solving function, called as obj.solve();

scvx.print_intro(obj.name,obj.ctrl.algo) 

scaling = set_scales(obj);

for iter = 1:obj.ctrl.iter_max
    
    % solve step
    varxu = solve_socp(obj,scaling,iter);
    
    % re-convexify along new solution (needed to get new defects)
    obj.convexify();
           
    % perform update step
    if (strcmp(obj.ctrl.algo,'ptr'))
        reject = false;
        change = '';
    else
        [reject,change] = obj.update_step();
    end
        
    % check convergence and exit if done
    converged = check_convergence(obj,varxu);
    
    % print some stuff
    print_status(obj,iter,varxu,reject,change);
    if (converged)
        break;
    end  
end

if (nargout>0)
    varargout{1} = get_exitcode(obj,iter,converged);
else
    get_exitcode(obj,iter,converged);
end

end

