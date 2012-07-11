function [] = basic(Tdir,infile,basename,tt)
%
% plots Salish simulations, a generic plot of the full domain
% 3/30/2011  Parker MacCready

% get file info
[G,S,T]=Z_get_basic_info(infile);

cmap = flipud(cbrewer('div','RdBu',32));

% plot salinity
s0 = nc_varget(infile,'salt',[0 S.N-1 0 0],[1 1 -1 -1]);
subplot(121)
pcolorcen(G.lon_rho,G.lat_rho,s0);
shading flat
cvec = [27 32]; caxis(cvec);
colorbar('Eastoutside');
% fix scaling
Z_dar;
% add labels
title('(a) Surface Salinity','fontweight','bold')
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
% add file info
Z_info(basename,tt,T,'lr'); 
% add coastline
Z_addcoast('detailed',Tdir.coast);
colormap(cmap);

% plot temperature
t0 = nc_varget(infile,'temp',[0 S.N-1 0 0],[1 1 -1 -1]);
subplot(122)
pcolorcen(G.lon_rho,G.lat_rho,t0);
shading flat
cvec = [10 13]; caxis(cvec);
colorbar('Eastoutside');
% fix scaling
Z_dar;
% add labels
title('(b) Surface Temperature','fontweight','bold')
xlabel('Longitude (deg)')
% add coastline
Z_addcoast('detailed',Tdir.coast);
% and velocity vectors
Z_velvec(infile,G,S,'lr')
colormap(cmap);
