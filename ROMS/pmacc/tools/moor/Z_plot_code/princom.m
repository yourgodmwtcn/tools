function []=princom(M)

% princom.m 2/5/2013 Parker MacCready
%
% plots the results of a mooring extraction

td = M.td;

u = Z_godin(M.ubar');
v = Z_godin(M.vbar');
mask = ~isnan(u);
u = u(mask); v = v(mask);
X = [u, v];

[COEFF,XX] = princomp(X);
ur = XX(:,1);
vr = XX(:,2);

figure; set(gcf,'position',[20 20 1250 750]); Z_fig;

subplot(121)
lon = M.run.grid.lon;
lat = M.run.grid.lat;
H = M.run.grid.H;
mask = M.run.grid.mask;
H(mask==0) = NaN;
pcolor(lon,lat,H);
shading interp;
Z_dar;
hold on
ph = plot(M.lon_rho,M.lat_rho,'pk','markersize',16);
set(ph,'markerfacecolor','y')
xlabel('Longitude')
ylabel('Latitude')
title([strrep(M.basename,'_',' '),' ',M.mloc],'fontweight','bold')

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



