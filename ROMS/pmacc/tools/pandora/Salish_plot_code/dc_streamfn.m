function [] = dc_streamfn(Tdir,infile,basename,tt)
% plots streamfunction - Deepak Cherian
% plots Salish simulations, a generic plot of the full domain
% 3/30/2011  Parker MacCready

% get file info
[G,S,T]=Z_get_basic_info(infile);

cmap = cbrewer('div','RdYlGn',32);
%[cmap] = cptcmap('saga-02','mapping','scaled');

% get vars
ubar      = double(squeeze(ncread(infile,'ubar',[1 1 1],[Inf Inf Inf])));
vbar      = double(squeeze(ncread(infile,'vbar',[1 1 1],[Inf Inf Inf])));
zeta      = nc_varget(infile,'zeta',[0 0 0],[Inf Inf Inf]);

% set up grid structure
[grid.xu,grid.yu] = toxy(G.lon_u(1,:),G.lat_u(:,1));
[grid.xv,grid.yv] = toxy(G.lon_v(1,:),G.lat_v(:,1));
[grid.xr,grid.yr] = toxy(G.lon_rho(1,:),G.lat_rho(:,1));

% calculate streamfunction
[psi,xpsi,ypsi] = streamfn_cgrid(grid,ubar,vbar,(zeta+G.h)');

[xpsi,ypsi] = tolatlon(xpsi,ypsi);

% plot
%subplot(121)
imagescnan(xpsi,ypsi,psi'); shading flat; hold on
colorbar('Southoutside');
% fix scaling
Z_dar;
% add labels
title('Transport Streamfunction','fontweight','bold')
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
Z_addcoast('detailed',Tdir.coast);
Z_velvec(infile,G,S,'mr')
% add file info
Z_info(basename,tt,T,'lr'); 
% add coastline

colormap(cmap);

function [x,y] = toxy(lon,lat)

    % Change bathy long,lat to x,y.
    latref = 48.4610;
    lonref = -122.955;

    latscale = 111190.29;
    lonscale = 74625.33;

    x = (lon-lonref)*lonscale/1000; % x-coord in km
    y = (lat-latref)*latscale/1000; % y-coord in km
    
function [lon,lat] = tolatlon(x,y)

    % Change bathy long,lat to x,y.
    latref = 48.4610;
    lonref = -122.955;

    latscale = 111190.29;
    lonscale = 74625.33;
    
    lon = lonref + 1000*x/lonscale;
    lat = latref + 1000*y/latscale;
