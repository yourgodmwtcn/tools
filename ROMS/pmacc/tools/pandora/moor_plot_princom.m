% moor_plot_princom.m 2/15/2012 Parker MacCready
%
% plots the results of a mooring extraction,
% finding the pricipal component direction of ubar, vbar

clear; close all;
[Tdir] = pan_start;
indir = [Tdir.pan_results,'moorings/'];
[fn,pth]=uigetfile([indir,'*.mat'], ...
    'Select mooring file...');
load([pth,fn]);
td = mod_moor.run.his.nctime;

u = mod_moor.ubar';
v = mod_moor.vbar';
X = [u, v];

[COEFF,XX] = princomp(X);
ur = XX(:,1);
vr = XX(:,2);

figure; set(gcf,'position',[20 20 1250 750]); Z_fig;

subplot(121)
lon = mod_moor.run.grid.lon;
lat = mod_moor.run.grid.lat;
H = mod_moor.run.grid.H;
mask = mod_moor.run.grid.mask;
H(mask==0) = NaN;
pcolor(lon,lat,H);
shading interp;
Z_dar;
Z_addcoast('combined',Tdir.coast);
hold on
ph = plot(mod_moor.p_lon,mod_moor.p_lat,'pk','markersize',16);
set(ph,'markerfacecolor','y')
xlabel('Longitude')
ylabel('Latitude')
title(strrep(fn,'_',' '))

vs = 1.2*max(sqrt(u.^2 + v.^2));

subplot(222)
plot(u,v,'*r')
axis equal
axis(vs*[-1 1 -1 1]);
grid on
hold on
c1 = COEFF(1,1); c2 = COEFF(2,1);
% theta is the angle (ccw positive) to rotate
% the axes to get to PC coordinates
theta = 180*atan2(c2,c1)/pi;
plot(vs*[-COEFF(1,1) COEFF(1,1)], ...
    vs*[-COEFF(2,1) COEFF(2,1)],'-k','linewidth',2);
title(['\theta = ',num2str(theta),' degrees'])
xlabel('U (m s^{-1})')
ylabel('V (m s^{-1})')

subplot(224)
plot(ur,vr,'*b')
axis equal
axis(vs*[-1 1 -1 1]);
grid on
xlabel('U_{ROTATED} (m s^{-1})')
ylabel('V_{ROTATED} (m s^{-1})')
title('Coordinates rotated ccw by \theta')



