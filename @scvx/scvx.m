classdef scvx < matlab.mixin.Copyable
    %scvx Basic implementation of Successive Convexification
    %
    % T. Reynolds -- RAIN Lab
        
    % Properties set by the user
    properties (SetAccess = public)
        name
        nx
        nu
        np
        auxdata
        bnds
        ctrl
        iguess
    end
    
    % Dynamics and Linearization function
    properties (SetAccess = public, Hidden = true)
        cost
        dynamics
        linearize
    end
    
    % Observable properties that cannot be set by the user
    properties (SetAccess = private, GetAccess = public)
        output
    end
    
    % Hidden properties that are set internally but still available to the
    % user.
    properties (SetAccess = private, Hidden = true)
        final_time_free
        initial_time_free
        best_output
        defects
        last_J
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%                        PUBLIC METHODS                       %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access=public)
        function obj = scvx(prob_name,nx,nu,np)
            %SCVX default constructor
            %   First input is the problem name as a string. Inputs two,
            %   three and four are the number of states, controls and
            %   parameters, respectively.
            if (nargin>0)
                if (~ischar(prob_name))
                    error('Input prob_name must be a string.')
                end
                obj.name = prob_name;
                obj.nx = nx;
                obj.nu = nu;
                obj.np = np;
            else
                obj.name = 'default';
                obj.nx = [];
                obj.nu = [];
                obj.np = [];
            end

        end
        
        function init(obj,p,varargin)
            %INIT Computes the initial guess
            %   obj = obj.init(P) does the default interpolation with
            %   minimum control.
            %
            %   obj = obj.init(P,X) uses the user defined matrix X as the
            %   initial state guess
            %
            %   obj = obj.init(P,X,U) uses the user defined matrices X and U
            %   as the initial state/control guesses. If you only want to
            %   guess the control, simple use obj.init([],U)
            if (nargin<1)
                error('User must provide a guess for the parameters')
            end
            if (numel(p)~=obj.np)
                error('Size of initial guess for parameters is not right')
            end
            obj.auxdata.tau = linspace(0,1,obj.ctrl.N);
            obj.iguess.p = p;
            if (nargin<3)
                % initial state guess
                obj.iguess.x = iguess_x(obj);
                % initial control guess
                obj.iguess.u = iguess_u(obj);
            elseif (nargin==3)
                obj.iguess.x = varargin{1};
                obj.iguess.u = iguess_u(obj);
            elseif (nargin==4)
                if (isempty(varargin{1}))
                    obj.iguess.x = iguess_x(obj);
                    obj.iguess.u = varargin{2};
                else
                    obj.iguess.x = varargin{1};
                    obj.iguess.u = varargin{2};
                end
            end
            % initialize free final & initial time flags
            tf_min = obj.bnds.trgt.t_min;
            tf_max = obj.bnds.trgt.t_max;
            t0_min = obj.bnds.init.t_min;
            t0_max = obj.bnds.init.t_max;
            tf_free = (tf_min~=tf_max);
            t0_free = (t0_min~=t0_max);  
            if (tf_free && t0_free)
                obj.bnds.p_min = [ tf_min; t0_min; obj.bnds.p_min ];
                obj.bnds.p_max = [ tf_max; t0_max; obj.bnds.p_max ];
            elseif (tf_free)
                obj.bnds.p_min = [ tf_min; obj.bnds.p_min ];
                obj.bnds.p_max = [ tf_max; obj.bnds.p_max ];
            elseif (t0_free)
                obj.bnds.p_min = [ t0_min; obj.bnds.p_min ];
                obj.bnds.p_max = [ t0_max; obj.bnds.p_max ];
            end
            obj.final_time_free     = tf_free;
            obj.initial_time_free   = t0_free;
            
            % convexify along initial guess
            obj.convexify();
            
            % compute the initial (nonlinear) cost
            obj.last_J = obj.cost(obj.iguess.x(:),obj.iguess.u(:),...
                        obj.iguess.p(:),obj.ctrl.N) ...
                        + obj.ctrl.wvc * sum(obj.defects);
                        
            % initial trust region
            obj.output.tr  = 0.5;
            obj.output.trp = 0.1;
        end % init
        
        function attach(obj,x,u,p)
            %ATTACH Attaches a given solution 
            %
            % obj.attach(x,u,p) provides a method to place the solution 
            % (x,u,p) in the output struct, which is typically read-only.
            
            szex = numel(x);
            if (szex ~= (obj.nx*obj.ctrl.N))
                error('SCVX::ATTACH the number of elements in x must be nx times N');
            end
            szeu = numel(u);
            if (szeu ~= (obj.nu*obj.ctrl.N))
                error('SCVX::ATTACH the number of elements in u must be nu times N');
            end
            szep = numel(p);
            if (szep ~= (obj.np))
                error('SCVX::ATTACH the number of elements in p must be np');
            end
            obj.output.x = x(:);
            obj.output.u = u(:);
            obj.output.p = p(:); 
            
        end % attach(x,u,p)
        
        function attach_file(obj,filename)
            %ATTACH_FILE    Attaches a given solution from a file
            % 
            % obj.attach(filename) provides a method to place the solution 
            % data stored in file 'filename.txt' in the output struct, 
            % which is typically read-only. Calls obj.attach(x,u,p) after
            % reading the data from the file.
            
            try
                data = load(filename);
            catch
                error('SCVX::ATTACH could not load data from %s',filename);
            end
            
            % extract sizes and (x,u,p) from data
            sze = data(1:3);
            data = data(4:end);
            x = data(1:sze(1));
            u = data(sze(1)+(1:sze(2)));
            p = data(sze(1)+sze(2)+(1:sze(3)));

            obj.attach(x,u,p);
                 
        end % attach_file
        
        function scaling = get_scales(obj)
            %GET_SCALES     
            %
            % Returns the scaling matrices used in the solution by calling
            % the (internal) function set_scales and returning the output
            % as a structure
            
            scaling = set_scales(obj);
        end
                       
        % The solve method contains the main loop that runs SCvx to solve
        % the instance of the problem
        varargout = solve(obj);
        
        % Function that convexifies dynamics and constraints
        convexify(obj);
        discretize(obj);
        
        % Function to perform an open loop propagation of the current
        % control solution
        result = open_loop(obj);
        
        % function to save data to ROOT/data/ folder
        save_data(obj);
        
    end % methods (public)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%                        PRIVATE METHODS                      %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access=private)
        
        function x = iguess_x(obj)
            % compute initial state guess
            x_ic_min = obj.bnds.init.x_min;
            x_ic_max = obj.bnds.init.x_max;
            x_fc_min = obj.bnds.trgt.x_min;
            x_fc_max = obj.bnds.trgt.x_max;
            x_ic = zeros(numel(x_ic_min),1);
            x_fc = x_ic;
            x = zeros(numel(x_ic_min),obj.ctrl.N);
            for k = 1:numel(x_ic_min)
                % Take the mean of any min/max values that don't agree
                if (x_ic_min(k)~=x_ic_max(k))
                    x_ic(k) = 0.5 * (x_ic_min(k)+x_ic_max(k));
                else
                    x_ic(k) = x_ic_min(k);
                end
                % set any NaNs to be zero
                if (isnan(x_ic(k)))
                    x_ic(k) = 0.0;
                end
                % Repeat steps for final condition
                if (x_fc_min(k)~=x_fc_max(k))
                    x_fc(k) = 0.5 * (x_fc_min(k)+x_fc_max(k));
                else
                    x_fc(k) = x_fc_min(k);
                end
                if (isnan(x_fc(k)))
                    x_fc(k) = 0.0;
                end
                % Interpolate between x_ic(k) and x_fc(k), place result in
                % the kth row of x
                x(k,:) = linspace(x_ic(k),x_fc(k),obj.ctrl.N);
            end
        end %iguess_x
        
        function u = iguess_u(obj)
            % compute initial control guess using minimum values
            u_min = obj.bnds.u_min;
            u = zeros(numel(u_min),obj.ctrl.N);
            for k = 1:numel(u_min)
                u(k,:) = linspace(u_min(k),u_min(k),obj.ctrl.N);
            end
        end
        
        function set_best_output(obj)
            %SET_BEST_OUTPUT
            %
            % Places the current values in obj.output into the struct
            % obj.best_output to be recovered later if SCVX needs.
            fn = fieldnames(obj.output);
            for k = 1:numel(fn)
                if ( (~strcmp(fn{k},'data')) && (~strcmp(fn{k},'tr')) )
                    obj.best_output.(fn{k}) = obj.output.(fn{k});
                end
            end
        end
        
        function get_best_output(obj)
           %GET_BEST_OUTPUT
           %
           % Places the values in obj.best_output into the struct
           % obj.output to recover the last best iterate of SCVX.
           fn = fieldnames(obj.best_output);
           for k = 1:numel(fn)
               obj.output.(fn{k}) = obj.best_output.(fn{k});
           end
        end
              
        % What data to display as we proceed through the iterations
        function print_status(obj,iter,varxu,reject,change)
            %PRINT_STATUS
            %
            % print_stats(iter,varxu,reject,change) prints the results of
            % the most recent iterate to the command window.
      
            vc1 = norm(obj.output.vc,1);
            tr1 = norm(obj.output.tr,1);
            J   = obj.output.cost;
            feas = obj.output.feas;
            fprintf('Iter: %02d | VC: %2.2e | TR: %2.2e | cost: %+2.2e',...
                iter,vc1,tr1,J);
            fprintf(' | max_var: %2.2e',varxu)
            fprintf(' | feas: %d |  Reject: %d%s\n',feas,reject,change);
        end
        
        function converged = check_convergence(obj,varxu)
            %CHECK_CONVERGENCE
            %
            %   converged = check_convergence(obj,varxu) checks the value
            %   of varxu against the convergence tolerance given in
            %   obj.ctrl.cvrg_tol specified by the user
            if (varxu < obj.ctrl.cvrg_tol)
                converged = true;
            else
                converged = false;
            end
        end
        
        function varargout = get_exitcode(obj,iter,converged)
            %EXITCODE
            %
            % Possible exitcodes are:
            % SOLVER_FLAG - exitcode(1):
            %   0 : all iterations returned solved
            %   1 : solver optimal, but reported difficulties at least once
            %   2 : solver failed at least once and dual problem may be 
            %       infeasible.
            %   3 : at least one iteration appears to be primal infeasible
            %   4 : at least one unknown error in solver
            % SCVX_FLAG - exitcode(2)
            %   0 : converged and feasible
            %   1 : reached max iterations and IS feasible 
            %   2 : converged and IS NOT feasible
            %   3 : reached max iterations and IS NOT feasible
            
            itr_status = cell2mat(obj.output.data.status);
            feas = obj.output.feas;
            exitcode = NaN(2,1);
            if (any(itr_status))
                flags = find(itr_status);
                for k = 1:numel(flags)
                    exitcode(1) = itr_status(flags(k));
                    % break is an iteration failed/was infeasible. These
                    % take highest priority and will drive the returned
                    % value, even if downstream iterations are okay.
                    if ( (exitcode(1)==3) || (exitcode(1)==2) )
                        break;
                    end
                end
            else
                exitcode(1) = 0;
                if (feas)
                    if (converged)
                        exitcode(2) = 0;
                    elseif (iter==obj.ctrl.iter_max)
                        exitcode(2) = 1;
                    end
                else
                    if (converged)
                        exitcode(2) = 2;
                    elseif (iter==obj.ctrl.iter_max)
                        exitcode(2) = 3;
                    end
                end
            end
            fprintf('Done. Solver flag %d | SCVX flag %d\n',...
                        exitcode(1),exitcode(2))
            if (nargout>0)
                varargout{1} = exitcode;
            end
            obj.output.data.exitcode = exitcode;
        end
        
        % Function that solves the SOCP at each iteration
        varxu = solve_socp(obj,scaling,iter);
        
        % Compute scaling terms
        scaling = set_scales(obj);
        
        % SCVX update step
        [reject,change] = update_step(obj);
        
    end % methods (private)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%                        STATIC METHODS                       %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function s = saveobj(obj)
            s.name = obj.name;
            s.nx = obj.nx;
            s.nu = obj.nu;
            s.np = obj.np;
            s.auxdata = obj.auxdata;
            s.bnds = obj.bnds;
            s.ctrl = obj.ctrl;
            s.iguess = obj.iguess;
            s.output = obj.output;
            s.final_time_free = obj.final_time_free;
            s.initial_time_free = obj.initial_time_free;
            s.best_output = obj.best_output;
            s.defects = obj.defects;
            s.last_J = obj.last_J;
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%                    STATIC PRIVATE METHODS                   %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access=private, Static=true)
        
        function print_intro(name,algo)
           
            if (strcmp(algo,'ptr'))
                introstr = 'Solving %s with penalized trust region...\n';
            else
                introstr = 'Solving %s with SCvx...\n';
            end
            fprintf(introstr,name);
            
        end
        
        function X = rk4(fun,t,x_ic)
            %RK4
            %
            % An implementation of Runga-Kutta 4th order fixed time step
            % integration. Y = scvx.rk4(func,tspan,x0) where
            %   func  : a function handle in the style of Matlab's
            %           ode45 function, e.g., @(t,x)foo(t,x,varargin)
            %   tspan : the temporal points at which to evaluate the 
            %           integral 
            %   x0    : the initial condition
            %
            % T. Reynolds -- RAIN Lab
            
            
            N = length(t);
            n = length(x_ic);
            X = zeros(N,n);
            X(1,:) = x_ic';
            
            for k = 1:N-1
                tk  = t(k);
                xk  = X(k,:)';
                h   = t(k+1) - tk;
                
                k1 = feval(fun,tk,xk);
                k2 = feval(fun,tk+h/2,xk+h/2*k1);
                k3 = feval(fun,tk+h/2,xk+h/2*k2);
                k4 = feval(fun,tk+h,xk+h*k3);
                
                X(k+1,:) = (xk+h/6*(k1+2*k2+2*k3+k4))';
            end
        end
        
        function v = vectorize(A)
           %VECTORIZE   Vectorize a matrix 
           %
           % v = scvx.vectorize(A) vectorizes the matrix A by stacking its
           % columns. So if 
           %
           %    A = [ a1, a2, ..., aN ] is m by n
           %
           % then v = [ a1 ; a2 ; ... ; aN ] is m*n by 1.
           [m,n] = size(A);
           v = zeros(m*n,1);
           for k = 1:n
               v(m*(k-1)+1:m*k) = A(:,k);
           end
        end
        
        function A = set_block(A,B,i,j)
            %SET_BLOCK sets block (i,j) of matrix A to the value B
            %
            % Suppose B is an [n,m] matrix and A is a [p,q] matrix. It is
            % assumed that p/n and q/m are divisible such that
            % A is composed of exactly p/n and q/m blocks.
            %
            % Parameters
            % ----------
            % A : matrix
            %     The matrix for the block to be set in.
            % B : matrix
            %     The value to which the block is to be set.
            % i : int
            %     Row index.
            % j : int
            %     Column index.
            %
            % Returns
            % -------
            % A : matrix
            %     Same matrix as A (the parameter), with block (i,j)
            %     set to B.
            [n,m] = size(B);
            A((i-1)*n+1:i*n,(j-1)*m+1:j*m) = B;
        end
        
        function B = get_block(A,i,j,h,v)
            %GET_BLOCK gets block (i,j) of matrix A
            %
            % Get block (i,j) of matrix A assumed that it is composed
            % of [h,v] blocks.
            %
            % Parameters
            % ----------
            % A : matrix
            %     The matrix for the block to be set in.
            % i : int
            %     Row index.
            % j : int
            %     Column index.
            % h : int
            %     Number of blocks composing the rows of A.
            % j : int
            %     Number of blocks composing the columns of A.
            %
            % Returns
            % -------
            % B : matrix
            %    Value of block (i,j) of A.
            
            [p,q]   = size(A);
            n       = p/h;
            m       = q/v;
            B       = A((i-1)*n+1:i*n,(j-1)*m+1:j*m);
        end
        
    end % methods (private, static)
end % classdef scvx

