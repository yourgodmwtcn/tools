function [] = dc_ssh_vor(Tdir,infile,basename,tt)
%
% Modified Parker's basic.m code
% 09 July 2012 - Deepak Cherian

% get file info
[G,S,T]=Z_get_basic_info(infile);

cmap = flipud(cbrewer('div','RdYlGn',24));

% plot SSH
s0 = nc_varget(infile,'zeta',[0 0 0],[Inf Inf Inf]);
s0 = s0 - nanmean(s0(:));

h(1) = subplot(221)
pcolorcen(G.lon_rho,G.lat_rho,s0); shading flat;
cvec = [-0.35 0.35]; caxis(cvec);
hcbar = colorbar('Eastoutside');
hold on; contour(G.lon_rho,G.lat_rho,G.h,[50 100 150],'k');
% fix scaling
Z_dar;
% add labels
title('SSH (m)','fontweight','bold')
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
% add file info
Z_info(basename,tt,T,'lr'); 
% add coastline
Z_addcoast('detailed',Tdir.coast);
colormap(cmap);

% plot vorticity
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
    
h(2) = subplot(222)
pcolor(xvor,yvor,vor'); shading interp
cvec = [-5 5]*1E-3; caxis(cvec);
hcbar = colorbar('Eastoutside'); 
% fix scaling
Z_dar;
% add labels
title('Surface Vorticity (s^{-1})','fontweight','bold')
xlabel('Longitude (deg)')
% add coastline
Z_addcoast('detailed',Tdir.coast);
% and velocity vectors
Z_velvec(infile,G,S,'mr')
colormap(cmap);

% now do tides
h(3) = subplot(2,2,[3 4])

set(h(1),'Units','pixels')
linkprop(h,'Units');

for i=1:1
    set(h(i),'Position',get(h(i),'Position')+[0 -225 0 180]);
end
pos1 = get(h(3),'Position');
set(h(3),'Position',pos1 + [0 0 0 -175])

set(h(1),'Units','normalized')
linkprop(h,'Units');

plot_tides(T.time_datenum); 
beautify(fontSize);

