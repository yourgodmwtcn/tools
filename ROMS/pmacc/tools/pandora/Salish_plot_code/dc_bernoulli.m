function [] = dc_bernoulli(Tdir,infile,basename,tt)
%
% Modified Parker's basic.m code
% 10 July 2012 - Deepak Cherian

% get file info
[G,S,T] = Z_get_basic_info(infile);

cmap = flipud(cbrewer('div','RdYlGn',24));
cmap2 = flipud(autumn);flipud(lbmap(24,'RedBlue'));

%% plot gz
zt = nc_varget(infile,'zeta',[0 0 0],[Inf Inf Inf]);
zt = zt - nanmean(zt(:));
zt = 9.81*(zt - nanmin(zt(:)));

subplot(131)
pcolorcen(G.lon_rho,G.lat_rho,zt); %shading interp;
cvec = [0 8]; caxis(cvec);
hcbar = colorbar('SouthOutside');
hold on; contour(G.lon_rho,G.lat_rho,G.h,[50 100],'k');
% fix scaling
Z_dar;
% add labels
title('g\zeta','fontweight','bold')
ylabel('Latitude (deg)')
% add file info
Z_info(basename,tt,T,'lr'); 
% add coastline
Z_addcoast('detailed',Tdir.coast);
colormap(cmap2);
%xlabel('Longitude (deg)')
freezeColors; cbfreeze;
ax(1) = gca;

%% plot 1/2 (u^2+v^2) + g zeta
u = double((ncread(infile,'u',[1 1 S.N 1],[Inf Inf 1 Inf])));
v = double((ncread(infile,'v',[1 1 S.N 1],[Inf Inf 1 Inf])));
uavg = avgx(u.*u);
vavg = avgy(v.*v);
u2v2 = squeeze(0.5*(uavg(:,2:end-1) + vavg(2:end-1,:)));
subplot(132)
pcolorcen(G.lon_rho(2:end-1,2:end-1),G.lat_rho(2:end-1,2:end-1),(u2v2+zt(2:end-1,2:end-1)')'); %shading interp
caxis(cvec);
hcbar = colorbar('SouthOutside'); 
% fix scaling
Z_dar;
% add labels
title('g\zeta + 1/2 (u^2+v^2)','fontweight','bold')
% add coastline
Z_addcoast('detailed',Tdir.coast);
% and velocity vectors
%Z_velvec(infile,G,S,'lr')
colormap(cmap2);
%xlabel('Longitude (deg)')
freezeColors; cbfreeze;
ax(2) = gca;

%% plot vorticity
warning off
grid = roms_get_grid(infile,infile,0,1);
warning on
%t0 = nc_varget(infile,'temp',[0 S.N-1 0 0],[1 1 -1 -1]);
u      = double(squeeze(ncread(infile,'u',[1 1 S.N 1],[Inf Inf 1 Inf])));
v      = double(squeeze(ncread(infile,'v',[1 1 S.N 1],[Inf Inf 1 Inf])));

grid1.xv = grid.lon_v(1,:)';
grid1.yv = grid.lat_v(:,1);
grid1.zv = grid.z_v(end,1,1);

grid1.xu = grid.lon_u(1,:)';
grid1.yu = grid.lat_u(:,1);
grid1.zu = grid.z_u(end,1,1);

toUTM =  findstr(ncreadatt(infile,'lon_u','units'),'degree');

[vor,xvor,yvor,zvor] = vorticity_cgrid(grid1,u,v,toUTM);
    
subplot(133)
pcolor(xvor,yvor,vor'); shading interp
cvec = [-5 5]*1E-3; caxis(cvec);
hcbar = colorbar('SouthOutside'); 
% fix scaling
Z_dar;
% add labels
title('Surface Vorticity (s^{-1})','fontweight','bold')
%xlabel('Longitude (deg)')
% add coastline
Z_addcoast('detailed',Tdir.coast);
% and velocity vectors
Z_velvec(infile,G,S,'lr')
colormap(cmap);
%xlabel('Longitude (deg)')

ax(3) = gca;

linkaxes(ax);

function [um] = avgy(um)
    um = (um(:,1:end-1,:,:)+um(:,2:end,:,:))/2;

function [um] = avgx(um)
    um = (um(1:end-1,:,:,:)+um(2:end,:,:,:))/2;

function [um] = avgz(um)
    um = (um(:,:,1:end-1,:)+um(:,:,2:end,:))/2;
