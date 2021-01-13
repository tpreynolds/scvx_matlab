clear; close all
run('../utils/set_path.m')

% auxiliary parameters
auxdata         = struct;
auxdata.g       = -9.81;
auxdata.cd      = 0.2;
auxdata.Sd      = 0.02;
auxdata.rho     = 1.225;
auxdata.m       = 2.5;
auxdata.J       = 1;
auxdata.l       = 0.1;

% boundary conditions
bnds = struct;
bnds.init.t_min = 0.0;
bnds.init.t_max = 0.0;
bnds.init.x_min = [ 2.0; 15.0; 1.0; 0.0; 0.0; 0.0 ];
bnds.init.x_max = [ 2.0; 15.0; 1.0; 0.0; 0.0; 0.0 ];
bnds.trgt.t_min = 5.0;
bnds.trgt.t_max = 15.0;
bnds.trgt.x_min = [ 0.0; 2.0; 0.0; -0.0; 0.0; 0.0 ];
bnds.trgt.x_max = [ 0.0; 2.0; 0.0; -0.0; 0.0; 0.0 ];

bnds.x_min = [ -6.0; 0.0; -7.0; -7.0; -pi/2; -0.5 ];
bnds.x_max = [  6.0; 20.0; 7.0;  7.0;  pi/2;  0.5 ];
bnds.u_min = [ 10; 4; 4 ];
bnds.u_max = [ 50; 10; 10 ];
bnds.p_min = []; % excluding final time
bnds.p_max = []; % excluding final time
bnds.path.x     = {@(x)planarquad_constraint(x)}; 
bnds.path.x_cvx = {false}; 
bnds.path.u     = {};
bnds.path.u_cvx = {};

ctrl = struct;
ctrl.N          = 30;
ctrl.Nsub       = 10;
ctrl.iter_max   = 10;
ctrl.wvc        = 1e3;
ctrl.wtr        = 1e-2;
ctrl.wtrp       = 1e-3;
ctrl.cvrg_tol   = 1e-4;
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
O = scvx('planarquad',6,3,1);

% attach structs to the created object
O.cost      = @planarquad_cost;          % cost function
O.dynamics  = @planarquad_dynamics;      % main dynamics function
O.linearize = @planarquad_linearize;     % main linearization function
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
p0      = 5.0;

% call initializing function -- parameters come first, they are the only
% mandatory thing to guess.
O.init(p0);%,x0,u0);

% solve the problem or grab a solution from a file
exitcode = O.solve();

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
r   = X(1:2,:);
v   = X(3:4,:);
Tlim = [ 0 O.output.p ];

figure(1), clf
subplot(2,1,1), hold on, grid on, box on
plot(T,r,'ko','MarkerFaceColor','k','HandleVisibility','off') 
plot(result.t,result.x(1:2,:))
set(gca,'Xlim',Tlim)
ylabel('Position')
subplot(2,1,2), hold on, grid on, box on
plot(T,v,'ko','MarkerFaceColor','k','HandleVisibility','off') 
plot(result.t,result.x(3:4,:))
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

% define size/shape
H  = diag([1,1]);
% define center
xc = [1;8];

figure(4), clf, hold on, grid on, box on
plot(r(1,:),r(2,:),'ko','MarkerFaceColor','k')
plot(result.x(1,:),result.x(2,:),'b')
angles = linspace(0,2*pi);
circle = [ cos(angles); sin(angles) ];
obs    = sqrtm(H) * circle + xc;
fill(obs(1,:),obs(2,:),'r','FaceAlpha',0.5)
ylabel('$r_x$')
xlabel('$r_y$')
title('Trajectory')