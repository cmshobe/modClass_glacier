% simple 1d glacier model
% written june 12 2015 rsa and cole 

clear all
figure(1)
clf
figure(2)
clf
figure(3)
clf
figure(4)
clf
figure(5)
clf

%% initialize

% material properties
rho_i = 917;
g = 9.81;
A = 2.1e-16; % Pa-3 yr-1
slide_ratio = 0.05;% ratio of sliding speed to internal defm speed

% initial topography
% zbmax = 4000;
% zbmin = 2000;
% xstar = 7000;

dx = 100;
xmax = 30000;
x = dx/2:dx:xmax-(dx/2);
xedge = 0:dx:xmax;

%zb = zbmin+(zbmax-zbmin)*exp(-x/xstar);

% now a clear creek colorado case
zb0 = 3400;
slope0 = 0.0338;
slope0 = 0.03;
zb1 = 516;
xstar = 1687;
zbline = zb0 - (slope0*x);
zbexp = zb1 * exp(-x/xstar);
zb = zbline + zbexp;

% valley width as a function of distance downvalley
% turn this off by setting Wmin=W0;
W0 = 3000;
Wstar = 4000;
Wmin = 1700;
Wmin = W0; % disallows variable width
W = Wmin + (W0-Wmin)*exp(-x./Wstar);
Wedge = W(1:end-1)+0.5*diff(W); % interpolates valley width to cell edges
Wedge = [Wedge(1) Wedge Wedge(end)];

%
zbmin = min(zb);
zbmax = max(zb);

H = zeros(size(x));
%Hedge = zeros(size(xedge));
z = zb+H;

dt = 0.002;
tmax = 1000;

t = 0:dt:tmax;
imax = length(t);
nplots = 50;
tplot = tmax/nplots;

% now set up meteorology (some based upon brugger's suggestion in sawatch
% paper 2006)
ELA0 = 3300;
%ELA0=3400;
% allow sinusoidal variation in ela...turn off with sigma_ELA = 0
%sigma_ELA = 200;
sigma_ELA = 0;

%ELA = ELA0 + sigma_ELA*randn(size(t));
period = 300;
ELA = ELA0 + sigma_ELA*cos(2*pi*t/period);

dbdz = 0.01; % usually 0.01 m/y/m
bcap = 1.5; % m/yr (usually 2m/yr)
%bcap = 3; % if v v high disallows cap
b0 = dbdz*(z-ELA0);
b0 = min(b0,bcap);
minzb = find(zb==min(zb));
minb = b0(minzb);

nframe = 0;

%% run

for i = 1:imax
    
b = dbdz*(z-ELA(i)); % local net balance calculated at cell centers at ice surface
b = min(b,bcap);

Hedge = H(1:end-1)+0.5*diff(H); % interpolates ice thickness to cell edges
S = abs(diff(z)/dx); % slope of ice surface calculated at cell edges

Udef = (A/5).*((rho_i*g*S).^3).*(Hedge.^4); %mean deformation speed
Q = (A/5).*((rho_i*g*S).^3).*(Hedge.^5); % ice discharge due to internal deformation
Qsl = (slide_ratio * Udef) .* Hedge; % ice discharge attributable to sliding...v crude here
Q = Q + Qsl;
Q =[0 Q 0]; %takes care of boundary conditions

taub = rho_i*g*S.*Hedge; % basal shear stress

dHdt = b - (1./W).*(diff(Q.*Wedge)/dx); %continuity allowing width to vary
H = H + (dHdt*dt); %updates ice thickness
H = max(H,0);

z = zb+H; %updates topography
glacier = find(H>0);

term(i) = x(glacier(end)); % terminus position
V_ice(i) = sum(H(glacier).*W(glacier))*dx; % glacier volume
    
% now for some plotting
if rem(t(i),tplot)==0
    nframe = nframe + 1
    
    figure(1)
    %subplot(1,2,1)
    subplot('position',[0.1 0.1 0.7 0.85])
        %subplot('position',[left bottom width height])
    plot(x/1000,zb,'k','linewidth',2)
    hold on
    plot(x/1000,z)
    plot(x/1000,ELA0*ones(size(x)),'g--','linewidth',2)
    axis([0 xmax/1000 zbmin-50 zbmax+200])
    xlabel('Horizontal Distance (km)','fontname','arial','fontsize',18)
    ylabel('Elevation (m)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')
    %subplot(1,2,2)
    subplot('position',[0.85 0.1 0.1 0.85])
    plot(b,z,'r','linewidth',1)
    hold on
    plot(b0,zb,'k--','linewidth',2)
    plot(zeros(size(z)),z,'g--','linewidth',2)
    %axis([minb 1.5*bcap zbmin-50 zbmax+200])
    axis([minb 1.2*max(b) zbmin-50 zbmax+200])
    xlabel('b (m/yr)','fontname','arial','fontsize',18)
    %ylabel('Elevation (m)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')
    % now stamp times
    time=num2str(t(i));
    timetext=strcat(time,' years')
    text(minb+1,zbmax+100,timetext,'fontsize',14)
    hold off
    
    Qanal = cumsum(b)*dx; %analytic solution for ss ice discharge
    Qanal = max(Qanal,0);
    
    figure(2)
    plot(xedge/1000,Q/1000)
    hold on
    plot(x/1000,Qanal/1000,'g--')
    axis([0 xmax/1000 0 bcap*xmax/2000])
    xlabel('Horizontal Distance (km)','fontname','arial','fontsize',18)
    ylabel('Ice discharge (1000 m^2/yr)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')
    hold off
    pause(0.1)
    
    figure(4)
    subplot(2,1,1)
    plot(x/1000,H)
    hold on
    xlabel('Distance (km)','fontname','arial','fontsize',18)
    ylabel('Ice thickness (m)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')
    subplot(2,1,2)
    plot(xedge(2:end-1)/1000,taub/1e5)
    hold on
    xlabel('Distance (km)','fontname','arial','fontsize',18)
    ylabel('Taub (bars)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')
    
end

end



%% finalize

figure(3)
subplot(3,1,1)
plot(t,ELA,'r')
    xlabel('Time (years)','fontname','arial','fontsize',18)
    ylabel('ELA (m)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')
subplot(3,1,2)
plot(t,term/1000,'r')
    xlabel('Time (years)','fontname','arial','fontsize',18)
    ylabel('Terminus position (km)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')
subplot(3,1,3)
plot(t,V_ice/1e9,'b')
    xlabel('Time (years)','fontname','arial','fontsize',18)
    ylabel('Glacier volume (km^3)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')

%     figure(4)
% plot(x/1000,H,'k')
%     xlabel('Distance (km)','fontname','arial','fontsize',18)
%     ylabel('Ice thickness (m)','fontname','arial','fontsize',18)
%     set(gca,'fontsize',14,'fontname','arial')

taub_test = 0:1e4:3e5;
strain_rate = A * taub_test.^3;
figure(5)
plot(taub_test/1e5,strain_rate,'r','linewidth',2)
    xlabel('Taub (bars)','fontname','arial','fontsize',18)
    ylabel('Strain rate (/year)','fontname','arial','fontsize',18)
    set(gca,'fontsize',14,'fontname','arial')
