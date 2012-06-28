function [] = basic_sji(Tdir,infile,basename,tt)
%
% plots Salish simulations, a generic plot of the full domain
% 3/30/2011  Parker MacCready
%
% 6/1/2012 This one focuses on the San Juan Islands

% get file info
[G,S,T]=Z_get_basic_info(infile);

aa = [-123.8 -122.3 48.2 49.1];
cbl = 'Northoutside';
scvec = [28 32];
tcvec = [10 14];

% plot salinity
s0 = nc_varget(infile,'salt',[0 S.N-1 0 0],[1 1 -1 -1]);
subplot(121)
pcolor(G.lon_rho,G.lat_rho,s0);
shading interp
axis(aa);
caxis(scvec);
colorbar(cbl);
% fix scaling
Z_dar;
% add labels
title('(a) Surface Salinity','fontweight','bold')
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
% add file info
Z_info(basename,tt,T,'ur');
% add coastline
Z_addcoast('combined',Tdir.coast);

% plot temperature
t0 = nc_varget(infile,'temp',[0 S.N-1 0 0],[1 1 -1 -1]);
subplot(122)
pcolor(G.lon_rho,G.lat_rho,t0);
shading interp
axis(aa);
caxis(tcvec);
colorbar(cbl);
% fix scaling
Z_dar;
% add labels
title('(b) Surface Temperature','fontweight','bold')
xlabel('Longitude (deg)')
% add coastline
Z_addcoast('combined',Tdir.coast);
% and velocity vectors
Z_velvec(infile,G,S,'ur')
