clear; close all
run('../utils/set_path.m')

% auxiliary parameters
auxdata   = struct;
auxdata.id_a = 1:3;
auxdata.id_w = 4:6;
auxdata.J    = diag([1.2;0.8;0.2]);

% boundary conditions
bnds = struct;
bnds.init.t_min = 0.0;
bnds.init.t_max = 0.0;
bnds.init.x_min = [ 0.0; 0.0; 0.0; 0.0; 0.0; 0.0 ];
bnds.init.x_max = [ 0.0; 0.0; 0.0; 0.0; 0.0; 0.0 ];
bnds.trgt.t_min = 15.0;
bnds.trgt.t_max = 35.0;
bnds.trgt.x_min = [ pi/4; pi/8; pi/4; 0.0; 0.0; 0.0 ];
bnds.trgt.x_max = [ pi/4; pi/8; pi/4; 0.0; 0.0; 0.0 ];

bnds.x_min = [ -pi; -pi/2; -pi; -0.15; -0.15; -0.15 ];
bnds.x_max = [  pi;  pi/2;  pi;  0.15;  0.15;  0.15 ];
bnds.u_min = [ -8; -8; -8 ] * 1e-3;
bnds.u_max = [  8;  8;  8 ] * 1e-3;
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
ctrl.wvc        = 1e4;
ctrl.wtr        = 1e-1;
ctrl.wtrp       = 1e-2;
ctrl.cvrg_tol   = 1e-3;
ctrl.feas_tol   = 1e-2;
ctrl.rho0       = 0.0;
ctrl.rho1       = 0.1;
ctrl.rho2       = 0.9;
ctrl.tr_lb      = 0.001;
ctrl.tr_ub      = 10;
ctrl.alpha      = 2.0;
ctrl.beta       = 2.0;
ctrl.algo       = 'ptr';
ctrl.solver     = 'sdpt3';

% initialize scvx problem with a name and dimensions (name,nx,nu,np)
O = scvx('satellite',6,3,1);

% attach structs to the created object
O.cost      = @satellite_cost;          % cost function
O.dynamics  = @satellite_dynamics;      % main dynamics function
O.linearize = @satellite_linearize;     % main linearization function
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
p0      = 15.0;

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
a   = rad2deg( X(1:3,:) );
w   = rad2deg( X(4:6,:) );
Tlim = [ 0 O.output.p ];

figure(1), clf
subplot(2,1,1), hold on, grid on, box on
plot(T,a,'ko','MarkerFaceColor','k','HandleVisibility','off') 
plot(result.t,rad2deg( result.x(1:3,:)) )
set(gca,'Xlim',Tlim)
ylabel('Attitude [deg]')
subplot(2,1,2), hold on, grid on, box on
plot(T,w,'ko','MarkerFaceColor','k','HandleVisibility','off') 
plot(result.t,rad2deg( result.x(4:6,:)) )
set(gca,'Xlim',Tlim)
ylabel('Attitude Rates [deg/s]')

figure(2), clf, hold on, grid on, box on
plot(T,U(1,:),'-o','DisplayName','$u_x$')
plot(T,U(2,:),'-o','DisplayName','$u_y$')
plot(T,U(3,:),'-o','DisplayName','$u_z$')
plot([0 T(end)],[O.bnds.u_min(1) O.bnds.u_min(1)],'r--',...
    'DisplayName','bound')
plot([0 T(end)],[O.bnds.u_max(1) O.bnds.u_max(1)],'r--',...
    'HandleVisibility','off')
set(gca,'Xlim',Tlim,'Ylim',[O.bnds.u_min(1),O.bnds.u_max(1)])
ylabel('Thrust')
legend('show')
xlabel('Time')
%%
figure(3), clf, hold on, grid on, box on
xI = [1.0;0.0;0.0];
yB = [0.0;1.0;0.0];
angles = linspace(0,2*pi,200);
[az,el] = cart2sph(xI(1),xI(2),xI(3));
az = abs(az);
r = tan(deg2rad(40));
fov_zone = zeros(3,numel(angles));
for k = 1:numel(angles)
    fov_zone(1,k) = r * sin(angles(k));
    fov_zone(2,k) = r * cos(angles(k));
    fov_zone(3,k) = sqrt( (1+r^2) - fov_zone(1,k)^2 - fov_zone(2,k)^2 );
end
R_sez2ecef = rot(rad2deg(az),'z')' * rot(90-rad2deg(el),'y')';
fov        = R_sez2ecef * fov_zone;
[fov_az, fov_el, fov_alt] = cart2sph(fov(1,:),fov(2,:),fov(3,:));

yB_azel = zeros(2,numel(result.t));
for k = 1:numel(result.t)
    C_I2B = my_eul2dcm(result.x(1:3,k));
    yB_   = C_I2B' * yB;
    [temp_az,temp_el,temp_alt] = cart2sph(yB_(1),yB_(2),yB_(3));
    yB_azel(:,k) = rad2deg([ abs(temp_az); temp_el ]);
end
% plot exclusion zone
fill(rad2deg(fov_az),rad2deg(fov_el),'r',...
        'FaceAlpha',0.25,'EdgeColor','r','LineWidth',2)
plot(yB_azel(1,:),yB_azel(2,:),'b','LineWidth',2)
set(gca,'XLim',[-180,180],'XTick',-180:45:180,...
        'YLim',[-90,90],'YTick',-90:45:90);
xlabel('Azimuth [deg]','FontSize',16)
ylabel('Elevation [deg]','FontSize',16)