% 1d glacier model
% written by RSA, modified by WHA for fluid earth lab on 13 oct 2015

clear all
close all

%% USER INPUTS
% Ice rheology and material properties
rho_i = 917; % ice density [kg m^-3]
A = 2.1e-16; % Pa-3 yr-1

% Topography
zbMax = 4000; % maximum bedrock elevation [m]
zbSlope = -0.03; % bedrock surface slope [m/m]


% Mass balance / climate
zEla = 3700; % elevation of the equilibirum line elevation [m]
dbdz = 0.01; % mass balance lapse rate [m y^-1 m^-1] (usually 0.01 m/y/m or so)
bCap = 1; % maximum mass balance [m y^-1]

% Model domain
% Space
dx = 100; % spatial step [m]
xMax = 30000; % end of spatial domain [m]

% Time
dt = 0.0025; % timestep [y] (this needs to be small [2.5e-3 works] for numerical stability)
tmax = 750; % end of temporal domain [y]

% constants
g = 9.81; % acceleration due to gravity [m s^-2]

% plotting

nplots = 50;
tplot = tmax/nplots;

%% INITIALIZING

% Initializing space and time vectors
x = dx/2:dx:xMax-(dx/2); % x-coordinates of node centers [m]
xedge = 0:dx:xMax; % x-coordinate of node edges [m]
t = 0:dt:tmax;

zb = zbMax + zbSlope*x; % bedrock elevation [m]
H = zeros(size(x)); % ice thickness
zs = zb+H; % ice surface elevation [m]
b = dbdz*(zs-zEla); % mass balance profile [m y^-1]

%% run

for i = 1:length(t)
    % Recalculate mass balance
    b = dbdz*(zs-zEla); % local net balance calculated at cell centers at ice surface
    b = min(b,bCap); % limit to prescribed maximum mass balance

    Hedge = H(1:end-1)+0.5*diff(H); % interpolates ice thickness to cell edges
    alpha = diff(zs)/dx; % slope of ice surface calculated at cell edges

    Q = (2*A/5).*((rho_i*g*sin(-alpha)).^3).*(Hedge.^5); % ice discharge due to internal deformation
    Q =[0 Q 0]; %takes care of boundary conditions - no flux in at top and no flux out at bottom

    dHdt = b - (diff(Q)/dx); % time rate of change of ice thickness [m y^-1]
    H = H + (dHdt*dt); % update ice thickness
    H = max(H,0); % can't have negative thickness

    zs = zb+H; %updates surface elevation



    % now for some plotting
    if rem(t(i),tplot)==0
        disp(['Time: ' num2str(t(i))])
        figure(1)
        plot(x/1000,zb,'k-','linewidth',2)
        axis([0 x(end)/1000 3000 5000])
        hold on
        plot(x/1000,zs,'b-','linewidth',2)
        xlabel('Down-glacier distance [km]','fontsize',18)
        ylabel('Elevation[m]','fontsize',18)
        legend('Bedrock','Ice surface','location','northeast')
        set(gca,'fontsize',14)
        pause(0.01)

    end

end




