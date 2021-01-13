% initialize example
run('../scvx_startup');

% auxiliary parameters
auxdata         = struct;
auxdata.g       = -1.0;
auxdata.cd      = 0.5;
auxdata.Sd      = 0.05;
auxdata.rho     = 1.25;
auxdata.m       = 1;

% boundary conditions
bnds = struct;
bnds.init.t_min = 0.0;
bnds.init.t_max = 0.0;
bnds.init.x_min = [ 5.0; 6.0; 24.0; -4.0; -2.0; 0.0 ];
bnds.init.x_max = [ 5.0; 6.0; 24.0; -4.0; -2.0; 0.0 ];
bnds.trgt.t_min = 4.0;
bnds.trgt.t_max = 12.0;
bnds.trgt.x_min = [ 0.0; 0.0; 0.0; 0.0; 0.0; 0.0 ];
bnds.trgt.x_max = [ 0.0; 0.0; 0.0; 0.0; 0.0; 0.0 ];

bnds.x_min = [ -10.0; -10.0; -1.0; -5.0; -5.0; -5.0 ];
bnds.x_max = [  10.0;  10.0; 25.0;  5.0;  5.0;  5.0 ];
bnds.u_min = [ -2; -2; -2 ];
bnds.u_max = [  2;  2;  2 ];
bnds.p_min = []; % excluding final time
bnds.p_max = []; % excluding final time
bnds.path.x     = {};
bnds.path.x_cvx = {};
bnds.path.u     = {};
bnds.path.u_cvx = {};

ctrl = struct;
ctrl.N          = 20;
ctrl.Nsub       = 10;
ctrl.iter_max   = 20;
ctrl.wvc        = 1e2;
ctrl.wtr        = 1e-2;
ctrl.wtrp       = 1e-3;
ctrl.cvrg_tol   = 1e-2;
ctrl.feas_tol   = 1e-2;
ctrl.rho0       = 0.0;
ctrl.rho1       = 0.1;
ctrl.rho2       = 0.9;
ctrl.tr_lb      = 0.001;
ctrl.tr_ub      = 10;
ctrl.alpha      = 2.0;
ctrl.beta       = 2.0;
ctrl.algo       = 'scvx';
ctrl.solver     = 'ecos';

% initialize scvx problem with a name and dimensions (name,nx,nu,np)
O = scvx('pointmass',6,3,1);

% attach structs to the created object
O.cost      = @pointmass_cost;          % cost function
O.dynamics  = @pointmass_dynamics;      % main dynamics function
O.linearize = @pointmass_linearize;     % main linearization function
O.auxdata   = auxdata;
O.bnds      = bnds;
O.ctrl      = ctrl;

% generate initial guess
% x0      = zeros(O.nx,O.ctrl.N);
% x0(1,:) = linspace(bnds.x_max(1),0.5*(bnds.trgt.x_max(1)+bnds.trgt.x_min(1)),ctrl.N);
% x0(2,:) = linspace(6,0,ctrl.N);
% x0(3,:) = linspace(24,1,ctrl.N);
% x0(4,:) = linspace(-4,0,ctrl.N);
% x0(5,:) = linspace(-2,-0.1,ctrl.N);
% u0      = zeros(O.nu,O.ctrl.N);
% u0(1,:) = [ bnds.u_max(1)*ones(1,4), bnds.u_min(1)*ones(1,2) bnds.u_max(1)*ones(1,ctrl.N-6) ];
p0      = 8.0;

% call initializing function -- parameters come first, they are the only
% mandatory thing to guess.
O.init(p0);%,x0,u0);

% solve the problem or grab a solution from a file
exitcode = O.solve();
% O.attach_file('../cpp_hp/data/optimal_xup.txt');

% perform open loop propagation on the last control
result = O.open_loop();

%% Plot the results
set(0,'defaulttextinterpreter','latex','defaultAxesFontSize',16,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex',...
    'defaultLineLineWidth',1.5,...
    'defaultLineMarkerSize',4)

% computed (discrete) solution
X   = reshape( O.output.x, O.nx, O.ctrl.N );
U   = reshape( O.output.u, O.nu, O.ctrl.N );
T   = O.output.p .* O.auxdata.tau;
r   = X(1:3,:);
v   = X(4:6,:);
Tlim = [ 0 O.output.p ];

figure(1), clf
subplot(2,1,1), hold on, grid on, box on
plot(T,r,'ko','MarkerFaceColor','k','HandleVisibility','off') 
plot(result.t,result.x(1:3,:))
set(gca,'Xlim',Tlim)
ylabel('Position')
subplot(2,1,2), hold on, grid on, box on
plot(T,v,'ko','MarkerFaceColor','k','HandleVisibility','off') 
plot(result.t,result.x(4:6,:))
set(gca,'Xlim',Tlim)
ylabel('Velocity')

figure(2), clf, hold on, grid on, box on
plot(T,U(1,:),'-o','DisplayName','$u_x$')
plot(T,U(2,:),'-o','DisplayName','$u_y$')
plot(T,U(3,:),'-o','DisplayName','$u_z$')
plot([0 T(end)],[O.bnds.u_min(1) O.bnds.u_min(1)],'r--',...
    'DisplayName','bound')
plot([0 T(end)],[O.bnds.u_max(1) O.bnds.u_max(1)],'r--',...
    'HandleVisibility','off')
set(gca,'Xlim',Tlim,'Ylim',[floor(O.bnds.u_min(1)),ceil(O.bnds.u_max(1))])
ylabel('Thrust')
legend('show')
xlabel('Time')

% figure(4), clf, hold on, grid on, box on
% plot(r(1,:),r(2,:),'ko','MarkerFaceColor','k')
% plot(result.x(2,:),result.x(3,:),'b')
% ylabel('$z_{\mathcal{I}}$')
% xlabel('$y_{\mathcal{I}}$')
% title('Landing Trajectory')