function []=profile(M)

% profile.m 2/5/2013 Parker MacCready
%
% plots the results of a mooring extraction

td = M.td;
ys = datestr(td(1),'yyyy');
yn = str2num(ys);
td0 = td - datenum(yn,1,1,0,0,0);

z = M.z_rho; eta = M.zeta; zb = M.z_w(1,:);
u = M.u; v = M.v;
t = M.temp; s = M.salt;
% extrapolate to top and bottom
zz = [zb; z; eta];
uu = [0*zb; u; u(end,:)];
vv = [0*zb; v; v(end,:)];
tt = [t(1,:); t; t(end,:)];
ss = [s(1,:); s; s(end,:)];
[NR,NC] = size(ss);

td0 = repmat(td0,NR,1);

figure; set(gcf,'position',[20 20 1400 900]); Z_fig(16);

subplot(411)
pcolor(td0,zz,uu)
shading interp
caxis([-0.5 0.5]);
colorbar('eastoutside')
ylabel('Z (m)')
axis([td0(1) td0(end) zz(1,1) 5]);
set(gca,'xticklabel',[]);
[xt,yt] = Z_lab('ll');
text(xt,yt,'U (m s^{-1}) ')
title([strrep(M.basename,'_',' '),' ',M.mloc],'fontweight','bold')

subplot(412)
pcolor(td0,zz,vv)
shading interp
caxis([-0.5 0.5]);
colorbar('eastoutside')
ylabel('Z (m)')
axis([td0(1) td0(end) zz(1,1) 5]);
set(gca,'xticklabel',[]);
[xt,yt] = Z_lab('ll');
text(xt,yt,'V (m s^{-1}) ')

subplot(413)
pcolor(td0,zz,tt)
shading interp
caxis([4 17]);
colorbar('eastoutside')
hold on
size(td0)
size(zz)
size(tt)

%contour(td0,zz,tt,[4:.5:25],'-k');
pcolor(td0,zz,tt); shading flat;
ylabel('Z (m)')
axis([td0(1) td0(end) zz(1,1) 5]);
set(gca,'xticklabel',[]);
[xt,yt] = Z_lab('ll');
text(xt,yt,'Temperature (\circC) ','color','w')

subplot(414)
pcolor(td0,zz,ss)
shading interp
caxis([28 33]);
colorbar('eastoutside')
hold on
%contour(td0,zz,ss,[0:.2:36],'-k');
ylabel('Z (m)')
axis([td0(1) td0(end) zz(1,1) 5]);
[xt,yt] = Z_lab('ll');
text(xt,yt,'Salinity ','color','w')
xlabel(['Yearday ',ys])

