clear; close all
run('../utils/set_path.m')

% auxiliary parameters
auxdata         = struct;
auxdata.g       = -1.0;
auxdata.alpha   = 1/(9.806*3);
auxdata.inertia = 0.5;

% boundary conditions
bnds = struct;
bnds.init.t_min = 0.0;
bnds.init.t_max = 0.0;
% bnds.init.x_min = [ 8984.0; 250.0; 433.0; -30.0; -15.0; -pi; 0.0 ];
% bnds.init.x_max = [ 8984.0; 250.0; 433.0; -30.0; -15.0; pi; 0.0 ];
% bnds.trgt.t_min = 20.0;
% bnds.trgt.t_max = 40.0;
% bnds.trgt.x_min = [ 8050.0; 0.0; 30.0; 0.0; -1; 0.0; 0.0 ];
% bnds.trgt.x_max = [ 8984.0; 0.0; 30.0; 0.0; -1; 0.0; 0.0 ];
bnds.init.x_min = [ 5.0; 6.0; 24.0; -4.0; -2.0; -pi; 0.0 ];
bnds.init.x_max = [ 5.0; 6.0; 24.0; -4.0; -2.0; pi; 0.0 ];
bnds.trgt.t_min = 4.0;
bnds.trgt.t_max = 12.0;
bnds.trgt.x_min = [ 2.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0 ];
bnds.trgt.x_max = [ 5.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0 ];

bnds.x_min = [ 2.0; -10.0; -1.0; -5.0; -5.0; -pi; -1.0 ];
bnds.x_max = [ 5.0; 10.0; 25.0; 5.0; 5.0; pi; 1.0 ];
bnds.u_min = [ 1.5; -0.1 ];
bnds.u_max = [ 6.5; 0.1 ];
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
ctrl.wtr        = 5e-2;
ctrl.wtrp       = 1e-2;
ctrl.cvrg_tol   = 1e-2;
ctrl.feas_tol   = 1e-2;
ctrl.rho0       = 0.0;
ctrl.rho1       = 0.1;
ctrl.rho2       = 0.9;
ctrl.tr_lb      = 0.001;
ctrl.tr_ub      = 10;
ctrl.alpha      = 2.0;
ctrl.beta       = 2.0;
ctrl.algo       = 'ptr';
ctrl.solver     = 'ecos';

% initialize scvx problem with a name and dimensions (name,nx,nu,np)
O = scvx('planarlanding',7,2,1);

% attach structs to the created object
O.cost      = @(x,u,p,N)( -x(7*(ctrl.N-1)+1) ); % cost function
O.dynamics  = @planarlanding_dynamics;        % main dynamics function
O.linearize = @planarlanding_linearize;       % main linearization function
O.auxdata   = auxdata;
O.bnds      = bnds;
O.ctrl      = ctrl;

% generate initial guess
x0      = zeros(O.nx,O.ctrl.N);
x0(1,:) = linspace(bnds.x_max(1),0.5*(bnds.trgt.x_max(1)+bnds.trgt.x_min(1)),ctrl.N);
x0(2,:) = linspace(6,0,ctrl.N);
x0(3,:) = linspace(24,1,ctrl.N);
x0(4,:) = linspace(-4,0,ctrl.N);
x0(5,:) = linspace(-2,-0.1,ctrl.N);
u0      = zeros(O.nu,O.ctrl.N);
u0(1,:) = [ bnds.u_max(1)*ones(1,4), bnds.u_min(1)*ones(1,2) bnds.u_max(1)*ones(1,ctrl.N-6) ];
p0      = 8.0;

% call initializing function -- parameters come first, they are the only
% mandatory thing to guess.
O.init(p0,x0,u0);

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
m   = X(1,:);
r   = X(2:3,:);
v   = X(4:5,:);
th  = X(6,:);
w   = X(7,:);
Tlim = [ 0 O.output.p ];

figure(1), clf
subplot(3,1,1), hold on, grid on, box on
plot(T,m,'ko','MarkerFaceColor','k','HandleVisibility','off')
plot(result.t,result.x(1,:))
set(gca,'Xlim',Tlim)
ylabel('Mass')
subplot(3,1,2), hold on, grid on, box on
plot(T,r,'ko','MarkerFaceColor','k','HandleVisibility','off') 
plot(result.t,result.x(2:3,:))
set(gca,'Xlim',Tlim)
ylabel('Position')
subplot(3,1,3), hold on, grid on, box on
plot(T,v,'ko','MarkerFaceColor','k','HandleVisibility','off') 
plot(result.t,result.x(4:5,:))
set(gca,'Xlim',Tlim)
ylabel('Velocity')

figure(2), clf
subplot(2,1,1), hold on, grid on, box on
plot(T,th,'ko','MarkerFaceColor','k','HandleVisibility','off')
plot(result.t,result.x(6,:))
set(gca,'Xlim',Tlim)
ylabel('Attitude angle')
subplot(2,1,2), hold on, grid on, box on
plot(T,w,'ko','MarkerFaceColor','k','HandleVisibility','off') 
plot(result.t,result.x(7,:))
set(gca,'Xlim',Tlim)
ylabel('Angular rate')

figure(3), clf
subplot(2,1,1), hold on, grid on, box on
plot(T,U(1,:),'b-o')
plot([0 T(end)],[O.bnds.u_min(1) O.bnds.u_min(1)],'r--')
plot([0 T(end)],[O.bnds.u_max(1) O.bnds.u_max(1)],'r--')
set(gca,'Xlim',Tlim,'Ylim',[0,ceil(O.bnds.u_max(1))])
ylabel('Thrust')
subplot(2,1,2), hold on, grid on, box on
plot(T,U(2,:),'b-o')
plot([0 T(end)],[O.bnds.u_min(2) O.bnds.u_min(2)],'r--')
plot([0 T(end)],[O.bnds.u_max(2) O.bnds.u_max(2)],'r--')
set(gca,'Xlim',Tlim,'Ylim',[1.5*(O.bnds.u_min(2)),1.5*(O.bnds.u_max(2))])
ylabel('Torque')
xlabel('Time')

figure(4), clf, hold on, grid on, box on
plot(r(1,:),r(2,:),'ko','MarkerFaceColor','k')
plot(result.x(2,:),result.x(3,:),'b')
ylabel('$z_{\mathcal{I}}$')
xlabel('$y_{\mathcal{I}}$')
title('Landing Trajectory')