% moor_plot_vs_depth.m 3/28/2012 Parker MacCready
%
% plots the results of a mooring extraction, focused on time series vs.
% depth

clear; %close all;
[Tdir] = pan_start;
indir = [Tdir.pan_results,'moorings/'];
[fn,pth]=uigetfile([indir,'*.mat'], ...
    'Select mooring file...');
load([pth,fn]);
td = mod_moor.run.his.nctime;
ys = datestr(td(1),'yyyy');
yn = str2num(ys);
td0 = td - datenum(yn,1,1,0,0,0);

z = mod_moor.z_rho; eta = mod_moor.zeta; zb = mod_moor.z_w(1,:);
u = mod_moor.u; v = mod_moor.v;
t = mod_moor.temp; s = mod_moor.salt;
% extrapolate to top and bottom
zz = [zb; z; eta];
uu = [0*zb; u; u(end,:)];
vv = [0*zb; v; v(end,:)];
tt = [t(1,:); t; t(end,:)];
ss = [s(1,:); s; s(end,:)];
[NR,NC] = size(ss);

td00 = ones(NR,1) * td0;


figure; set(gcf,'position',[20 20 1400 900]); Z_fig;

subplot(411)
pcolor(td0,zz,uu)
shading interp
colorbar('eastoutside')
ylabel('Z (m)')
axis([td0(1) td0(end) zz(1,1) 5]);
set(gca,'xticklabel',[]);
[xt,yt] = Z_lab('ll');
text(xt,yt,'U (m s^{-1}) ')
%
name = strrep(fn,'_','\_'); name = strrep(name,'.mat','');
title(name,'fontweight','bold')

subplot(412)
pcolor(td0,zz,vv)
shading interp
colorbar('eastoutside')
ylabel('Z (m)')
axis([td0(1) td0(end) zz(1,1) 5]);
set(gca,'xticklabel',[]);
[xt,yt] = Z_lab('ll');
text(xt,yt,'V (m s^{-1}) ')

subplot(413)
pcolor(td0,zz,tt)
shading interp
colorbar('eastoutside')
hold on
contour(td00,zz,tt,[4:.5:25],'-k');
ylabel('Z (m)')
axis([td0(1) td0(end) zz(1,1) 5]);
set(gca,'xticklabel',[]);
[xt,yt] = Z_lab('ll');
text(xt,yt,'Temperature (\circC) ')

subplot(414)
pcolor(td0,zz,ss)
shading interp
colorbar('eastoutside')
hold on
contour(td00,zz,ss,[0:.1:36],'-k');
ylabel('Z (m)')
axis([td0(1) td0(end) zz(1,1) 5]);
[xt,yt] = Z_lab('ll');
text(xt,yt,'Salinity ')
%
xlabel(['Yearday ',ys])

