function [] = dc_bernoulli(Tdir,infile,basename,tt)
%
% Modified Parker's basic.m code
% 10 July 2012 - Deepak Cherian

% get file info
[G,S,T] = Z_get_basic_info(infile);

cmap = flipud(cbrewer('div','RdYlGn',24));

% plot gz
zt = nc_varget(infile,'zeta',[0 0 0],[Inf Inf Inf]);
zt = zt - nanmean(zt(:));
zt = 9.81*(zt - nanmin(zt(:)));
subplot(121)
pcolor(G.lon_rho,G.lat_rho,zt); shading interp;
cvec = [0 10]; caxis(cvec);
hcbar = colorbar('Eastoutside');
hold on; contour(G.lon_rho,G.lat_rho,G.h,[50 100 150],'k');
% fix scaling
Z_dar;
% add labels
title('g\zeta','fontweight','bold')
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
% add file info
Z_info(basename,tt,T,'lr'); 
% add coastline
Z_addcoast('detailed',Tdir.coast);
colormap(cmap);

% plot 1/2 u^2
u = double((ncread(infile,'u',[1 1 S.N 1],[Inf Inf 1 Inf])));
v = double((ncread(infile,'v',[1 1 S.N 1],[Inf Inf 1 Inf])));
uavg = avgx(u.*u);
vavg = avgy(v.*v);
u2v2 = squeeze(0.5*(uavg(:,2:end-1) + vavg(2:end-1,:)));
subplot(122)
pcolor(G.lon_rho(2:end-1,2:end-1),G.lat_rho(2:end-1,2:end-1),u2v2'); shading interp
cvec = [0 10]; caxis(cvec);
hcbar = colorbar('Eastoutside'); 
% fix scaling
Z_dar;
% add labels
title('1/2 (u^2+v^2)','fontweight','bold')
xlabel('Longitude (deg)')
% add coastline
Z_addcoast('detailed',Tdir.coast);
% and velocity vectors
Z_velvec(infile,G,S,'lr')
colormap(cmap);

function [um] = avgy(um)
    um = (um(:,1:end-1,:,:)+um(:,2:end,:,:))/2;

function [um] = avgx(um)
    um = (um(1:end-1,:,:,:)+um(2:end,:,:,:))/2;

function [um] = avgz(um)
    um = (um(:,:,1:end-1,:)+um(:,:,2:end,:))/2;
