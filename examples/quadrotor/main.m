% initialize example
run('../scvx_startup');

% auxiliary parameters
auxdata           = struct;
auxdata.g         = -9.81;
auxdata.u_nrm_min = 0.6;
auxdata.u_nrm_max = 23.2;
auxdata.tilt_max  = deg2rad(60);
auxdata.id_r = 1:3;
auxdata.id_v = 4:6;
auxdata.id_w = 1:3;
auxdata.id_G = 4;

% boundary conditions
bnds = struct;
bnds.init.t_min = 0.0;
bnds.init.t_max = 0.0;
bnds.init.x_min = [ 0.0; 0.0; 0.0; 0.0; 0.0; 0.0 ];
bnds.init.x_max = [ 0.0; 0.0; 0.0; 0.0; 0.0; 0.0 ];
bnds.trgt.t_min = 2.0;
bnds.trgt.t_max = 3.0;
bnds.trgt.x_min = [ 2.5; 6.0; 0.0; 0.0; 0.0; 0.0 ];
bnds.trgt.x_max = [ 2.5; 6.0; 0.0; 0.0; 0.0; 0.0 ];

bnds.x_min = [ -10.0; -10.0; -10.0; -8.0; -8.0; -8.0 ];
bnds.x_max = [ 10.0; 10.0; 10.0; 8.0; 8.0; 8.0 ];
bnds.u_min = [ -23.2; -23.2; 0.0; 0.0 ];
bnds.u_max = [ 23.2; 23.2; 23.2; 23.2 ];
bnds.p_min = []; % excluding final time
bnds.p_max = []; % excluding final time
bnds.path.x     = {}; 
bnds.path.x_cvx = {}; 
bnds.path.u     = {@(u)quadrotor_control_constraints(u,auxdata)};
bnds.path.u_cvx = {true};

ctrl = struct;
ctrl.N          = 20;
ctrl.Nsub       = 10;
ctrl.iter_max   = 20;
ctrl.wvc        = 1e2;
ctrl.wtr        = 1e-2;
ctrl.wtrp       = 1e-3;
ctrl.cvrg_tol   = 1e-3;
ctrl.feas_tol   = 1e-2;
ctrl.rho0       = 0.0;
ctrl.rho1       = 0.1;
ctrl.rho2       = 0.7;
ctrl.tr_lb      = 0.001;
ctrl.tr_ub      = 10;
ctrl.alpha      = 2.0;
ctrl.beta       = 2.0;
ctrl.algo       = 'scvx';
ctrl.solver     = 'ecos';

% initialize scvx problem with a name and dimensions (name,nx,nu,np)
O = scvx('quadrotor',6,4,1);

% attach structs to the created object
O.cost      = @quadrotor_cost;          % cost function
O.dynamics  = @quadrotor_dynamics;      % main dynamics function
O.linearize = @quadrotor_linearize;     % main linearization function
O.auxdata   = auxdata;
O.bnds      = bnds;
O.ctrl      = ctrl;

% generate initial guess
p0      = 2.0;

% call initializing function -- parameters come first, they are the only
% mandatory thing to guess. If no x0, u0 are passed, they will be
% automatically guessed
O.init(p0);

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
u   = U(1:3,:);
G   = U(4,:);
u_nrm = zeros(size(G));
u_ang = zeros(size(G));
for k = 1:O.ctrl.N
    u_nrm(k) = norm(u(:,k));
    u_ang(k) = acosd(u(3,k)./u_nrm(k));
end
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

figure(2), clf
subplot(2,2,[1,2]), hold on, grid on, box on
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
subplot(2,2,3), hold on, grid on, box on
plot(T,G,'r--')
plot(T,u_nrm,'k')
ylabel('Thrust Magnitude')
xlabel('Time')
subplot(2,2,4), hold on, grid on, box on
plot(T,u_ang,'k')
plot([0 T(end)],rad2deg([O.auxdata.tilt_max(1) O.auxdata.tilt_max(1)]),'r--',...
    'DisplayName','bound')
plot([0 T(end)],-rad2deg([O.auxdata.tilt_max(1) O.auxdata.tilt_max(1)]),'r--',...
    'HandleVisibility','off')
ylabel('Tilt Angle')
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