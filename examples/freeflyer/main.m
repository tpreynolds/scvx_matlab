% initialize example
run('../scvx_startup');

% auxiliary parameters
auxdata           = struct;
auxdata.r_nrm_max = 12.0;
auxdata.v_nrm_max = 0.4;
auxdata.w_nrm_max = 1.0;
auxdata.T_nrm_max = 72e-3;
auxdata.M_nrm_max = 2e-3;
auxdata.id_r = 1:3;
auxdata.id_v = 4:6;
auxdata.id_q = 7:10;
auxdata.id_w = 11:13;
auxdata.id_T = 1:3;
auxdata.id_M = 4:6;
auxdata.mass = 7.2;
auxdata.inertia = diag([0.1083,0.1083,0.1083]);

% boundary conditions
bnds = struct;
bnds.init.t_min = 0.0;
bnds.init.t_max = 0.0;
bnds.init.x_min = [ 7.1; 0.4; 5.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0 ];
bnds.init.x_max = [ 7.3; 0.5; 5.5; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0 ];
bnds.trgt.t_min = 90.0;
bnds.trgt.t_max = 100.0;
bnds.trgt.x_min = [ 11.3; 6.0; 4.5; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0 ];
bnds.trgt.x_max = [ 11.3; 6.0; 4.5; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 1.0; 0.0; 0.0; 0.0 ];

bnds.x_min = [ -15.0; -15.0; -15.0; -0.5; -0.5; -0.5; -1.0; -1.0; -1.0; -1.0; -1.0; -1.0; -1.0 ];
bnds.x_max = [  15.0;  15.0;  15.0;  0.5;  0.5;  0.5;  1.0;  1.0;  1.0;  1.0;  1.0;  1.0;  1.0 ];
bnds.u_min = [ -8e-2; -8e-2; -8e-2; -1e-2; -1e-2; -1e-2 ];
bnds.u_max = [  8e-2;  8e-2;  8e-2;  1e-2;  1e-2;  1e-2 ];
bnds.p_min = []; % excluding final time
bnds.p_max = []; % excluding final time
bnds.path.x     = {}; 
bnds.path.x_cvx = {}; 
bnds.path.u     = {@(u)freeflyer_control_constraints(u,auxdata)};
bnds.path.u_cvx = {true};

ctrl = struct;
ctrl.N          = 50;
ctrl.Nsub       = 15;
ctrl.iter_max   = 15;
ctrl.wvc        = 1e4;
ctrl.wtr        = 1e-1;
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
ctrl.solver     = 'sdpt3';

% initialize scvx problem with a name and dimensions (name,nx,nu,np)
O = scvx('freeflyer',13,6,1);

% attach structs to the created object
O.cost      = @freeflyer_cost;          % cost function
O.dynamics  = @freeflyer_dynamics;      % main dynamics function
O.linearize = @freeflyer_linearize;     % main linearization function
O.auxdata   = auxdata;
O.bnds      = bnds;
O.ctrl      = ctrl;

% generate initial guess
p0      = 100.0;

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
r   = X(auxdata.id_r,:);
v   = X(auxdata.id_v,:);
F   = U(auxdata.id_T,:);
M   = U(auxdata.id_M,:);
F_nrm = zeros(O.ctrl.N,1);
M_nrm = zeros(O.ctrl.N,1);
for k = 1:O.ctrl.N
    F_nrm(k) = norm(F(:,k));
    M_nrm(k) = norm(M(:,k));
end
Tlim = [ 0 O.output.p ];

figure(1), clf
subplot(2,1,1), hold on, grid on, box on
plot(T,r,'ko','MarkerFaceColor','k','HandleVisibility','off') 
plot(result.t,result.x(auxdata.id_r,:))
set(gca,'Xlim',Tlim)
ylabel('Position')
subplot(2,1,2), hold on, grid on, box on
plot(T,v,'ko','MarkerFaceColor','k','HandleVisibility','off') 
plot(result.t,result.x(auxdata.id_v,:))
set(gca,'Xlim',Tlim)
ylabel('Velocity')

figure(2), clf
subplot(2,2,1), hold on, grid on, box on
plot(T,F(1,:),'-o','DisplayName','$T_x$')
plot(T,F(2,:),'-o','DisplayName','$T_y$')
plot(T,F(3,:),'-o','DisplayName','$T_z$')
plot(Tlim,[O.bnds.u_min(1) O.bnds.u_min(1)],'r--',...
    'DisplayName','bound')
plot(Tlim,[O.bnds.u_max(1) O.bnds.u_max(1)],'r--',...
    'HandleVisibility','off')
% set(gca,'Xlim',Tlim,'Ylim',[floor(O.bnds.u_min(1)),ceil(O.bnds.u_max(1))])
ylabel('Thrust')
legend('show')
xlabel('Time')
subplot(2,2,2), hold on, grid on, box on
plot(T,F_nrm,'k')
plot(Tlim,auxdata.T_nrm_max.*[1,1],'r--')
ylabel('Thrust Magnitude')
xlabel('Time')
%
subplot(2,2,3), hold on, grid on, box on
plot(T,M(1,:),'-o','DisplayName','$M_x$')
plot(T,M(2,:),'-o','DisplayName','$M_y$')
plot(T,M(3,:),'-o','DisplayName','$M_z$')
plot(Tlim,[O.bnds.u_min(4) O.bnds.u_min(4)],'r--',...
    'DisplayName','bound')
plot(Tlim,[O.bnds.u_max(4) O.bnds.u_max(4)],'r--',...
    'HandleVisibility','off')
ylabel('Torque')
xlabel('Time')
%
subplot(2,2,4), hold on, grid on, box on
plot(T,M_nrm,'k')
plot(Tlim,auxdata.M_nrm_max.*[1,1],'r--')
ylabel('Torque Magnitude')
xlabel('Time')

% define size/shape
H  = diag([1,1]);
% define center
xc = [1;8];

figure(4), clf, hold on, grid on, box on
plot3(r(1,:),r(2,:),r(3,:),'ko','MarkerFaceColor','k')
plot3(result.x(1,:),result.x(2,:),result.x(3,:),'b')
% angles = linspace(0,2*pi);
% circle = [ cos(angles); sin(angles) ];
% obs    = sqrtm(H) * circle + xc;
% fill(obs(1,:),obs(2,:),'r','FaceAlpha',0.5)
ylabel('$r_x$')
xlabel('$r_y$')
zlabel('$r_z$')
title('Trajectory')