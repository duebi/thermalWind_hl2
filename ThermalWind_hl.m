% ThermalWind_hl

load data/sg146m11data % load all seaglider data
load data/colorbrewer_anom

% 1. Isolate data from meridional transect
merid1 = dived.dive >= 50 & dived.dive <= 137 & dived.dive ~= 106; % indeces for first meridional transect
tw.lat = dived.lat(merid1); tw.lon = dived.lon(merid1); tw.date = dived.date(merid1);
tw.depth = sgd.depth; tw.sig = sgd.sig(:,merid1);

% 2. Compute horizontal distance
tw.dist = NaN*tw.lat;
tw.dist(1) = 0;
for i = 2:length(tw.lat)
    tw.dist(i) = tw.dist(i-1) + vdist(tw.lat(i-1),nanmean(tw.lon),tw.lat(i),nanmean(tw.lon)); % distance in meters
end

% 3. Reinterpolate data on regular grid (step < than maximum array step)
xstep = 500; % in meters
zstep = tw.depth(2) - tw.depth(1); % in meters
xint = 0:xstep:225000; % regular array of meridional distance 
zint = tw.depth;
latint = interp1(tw.dist,tw.lat,xint);
[xgrid,zgrid] = meshgrid(xint,zint);
sigint = interp2(tw.dist,tw.depth,tw.sig,xgrid,zgrid);

% 4. Smoothed density array (based on running mean)
sigint_smooth = NaN*sigint;
for i = 1:length(tw.depth)
    sigint_smooth(i,:) = nanrunmean(sigint(i,:)',31);
end
%{
%plot smoothed and raw data
contour(xint,zint,sigint,23:0.5:26,'edgecolor','k')
hold on, [c,h]=contour(xint,zint,sigint_smooth,23:0.5:26,'edgecolor','r'), hold off
clabel(c,h)
set(gca,'ydir','rev','Fontsize',16), ylim([0 350])
xlabel('meridional distance (m)'), ylabel('depth (m)')
%}

% 5.Compute vertical shear from horizontal density gradients (thermal wind approximation)
f = 2*7.2921e-5*sind(latint); %s-1
f = repmat(f,length(tw.depth),1);
rho_y = (sigint_smooth(:,3:end)-sigint_smooth(:,1:end-2))/(2*xstep);% kg m-4 (central finite difference)
rho_y = [NaN(length(tw.depth),1) rho_y NaN(length(tw.depth),1)];
u_z =   (9.8*rho_y)./((1000+sigint_smooth).*f); % shear of u

%{
% contour of meridional shear
contourf(xint,zint,u_z,-0.0305:0.001:0.0305,'edgecolor','none')
cb = colorbar; title(cb,'Shear of U (s-1)'); set(cb,'Fontsize',16)
caxis([-0.0055 0.0055])
hold on, contour(xint,zint,sigint_smooth,23:0.5:26,'edgecolor','k'), hold off
set(gca,'ydir','rev','Fontsize',16), ylim([0 350])
xlabel('meridional distance (m)'), ylabel('depth (m)')
colormap(anom_map2)
%}

% 6. Compute Brunt-Vaisala frequency
rho_z = (sigint_smooth(3:end,:)-sigint_smooth(1:end-2,:))/(2*zstep);
rho_z = [NaN(1,length(xint)); rho_z; NaN(1,length(xint))];
rho_z(rho_z<0) = NaN;
bvf = sqrt((9.8*rho_z)./(1000+sigint_smooth));
%{
contourf(xint,zint,bvf,0.00125:0.0025:0.04,'edgecolor','none')
cb = colorbar; title(cb,'Brunt-Vaisala freq. (s-1)'); set(cb,'Fontsize',16)
caxis([0 0.0225])
hold on, contour(xint,zint,sigint_smooth,23:0.5:26,'edgecolor','k'), hold off
set(gca,'ydir','rev','Fontsize',16), ylim([0 350])
xlabel('meridional distance (m)'), ylabel('depth (m)')
colormap(quant_blue) % 9 intervals
%}

% 6. Compute Richardson number
Ri = bvf./u_z;
