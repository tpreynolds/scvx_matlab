function convexify(obj)
%CONVEXIFY
%
% T. Reynolds -- RAIN Lab

% discretize along the current solution
obj.discretize();

% convexify and non-convex state constraints
nx_constraints_path = numel(obj.bnds.path.x);
for cnstr = 1:nx_constraints_path
   if (~obj.bnds.path.x_cvx{cnstr})
       for k = 1:obj.ctrl.N
           try
               xk = obj.output.x(obj.nx*(k-1)+(1:obj.nx));
           catch
               xk = obj.iguess.x(:,k);
           end
           [fk,Ak] = obj.bnds.path.x{cnstr}(xk);
           obj.bnds.path.x_lin{cnstr}{1}(k)   = fk - Ak*xk;
           obj.bnds.path.x_lin{cnstr}{2}(k,:) = Ak;
       end
   end
end

% TODO implement method to convexify any non-convex control constraints

end



