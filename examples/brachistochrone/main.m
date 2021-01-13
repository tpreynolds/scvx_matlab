% initialize example
run('../scvx_startup');

auxdata.g = 10;

bnds.init.t_min = 0.0;
bnds.init.t_max = 0.0;
bnds.init.x_min = [ 0.0; 0.0; 0.0 ];
bnds.init.x_max = [ 0.0; 0.0; 0.0 ];
bnds.trgt.t_min = 0.0;
bnds.trgt.t_max = 2.0;
bnds.trgt.x_min = [ 2.0; 2.0; -10.0 ];
bnds.trgt.x_max = [ 2.0; 2.0;  10.0 ];

bnds.x_min = [ 0.0; 0.0; -10.0 ];
bnds.x_max = [ 4.0; 4.0;  10.0 ];
bnds.u_min = [ -pi/2 ];
bnds.u_max = [ pi/2 ];
bnds.p_min = []; % excluding final time
bnds.p_max = []; % excluding final time
bnds.path.x     = {};
bnds.path.x_cvx = {};
bnds.path.u     = {};
bnds.path.u_cvx = {};

ctrl.N          = 30;
ctrl.Nsub       = 10;
ctrl.iter_max   = 25;
ctrl.wvc        = 1e3;
ctrl.wtr        = 1e-1;
ctrl.wtrp       = 1e-2;
ctrl.cvrg_tol   = 1e-2;
ctrl.feas_tol   = 1e-2;
ctrl.rho0       = 0.0;
ctrl.rho1       = 0.1;
ctrl.rho2       = 0.7;
ctrl.tr_lb      = 0.001;
ctrl.tr_ub      = 10;
ctrl.alpha      = 1.5;
ctrl.beta       = 2.0;
ctrl.algo       = 'ptr';
ctrl.solver     = 'sdpt3';

% initialize scvx problem with a name and dimensions (name,nx,nu,np)
O = scvx('brachistochrone',3,1,1);

% attach structs to the created object
O.cost      = @(x,u,p,N)( p(1) );
O.dynamics  = @brach_dynamics;    % main dynamics function
O.linearize = @brach_linearize;   % main linearization function
O.auxdata   = auxdata;
O.bnds      = bnds;
O.ctrl      = ctrl;

% generate initial guess
x0      = zeros(O.nx,O.ctrl.N);
x0(1,:) = linspace(0.0,2.0,ctrl.N);
x0(2,:) = linspace(0.0,2.0,ctrl.N);
x0(3,:) = linspace(0.0,0.0,ctrl.N);
u0      = linspace(0.0,pi/2,ctrl.N);
p0      = 2.0;

% call initializing function -- parameters come first, they are the only
% mandatory thing to guess.
O.init(p0,x0,u0);

% solve the problem
exitcode = O.solve();

% perform open loop propagation on the last control
result = O.open_loop();

%% Plot solution

% computed (discrete) solution
X   = reshape( O.output.x, O.nx, O.ctrl.N );
U   = reshape( O.output.u, O.nu, O.ctrl.N );
T   = O.output.p .* O.auxdata.tau;
xx   = X(1,:);
yy   = X(2,:);
vv   = X(3,:);
Tlim = [ 0 O.output.p ];

figure(1), clf, hold on, grid on, box on
pp = plot(result.t,result.x);
plot(T,X,'ko','MarkerSize',3,'HandleVisibility','off')
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$(x(t),y(t),v(t))$','Interpreter','LaTeX');
ll = legend('$x(t)$','$y(t)$','$v(t)$','Location','NorthWest');
set(pp,'LineWidth',1.25,'MarkerSize',8);
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',18,'Interpreter','LaTeX');
set(gca,'FontSize',16,'FontName','Times');

figure(2), hold on, grid on, box on
pp = plot(T,U,'k-o','MarkerSize',3);
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u(t)$','Interpreter','LaTeX');
ll = legend('$u(t)$','Location','NorthWest');
set(pp,'LineWidth',1.25,'MarkerSize',8);
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(ll,'FontSize',18,'Interpreter','LaTeX');
set(gca,'FontSize',16,'FontName','Times');