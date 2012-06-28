function [] = basic_zslice(Tdir,infile,basename,tt)
%
% plots Salish simulations, focused on the JdF eddy
% 11/6/2011  Parker MacCready
%
% uses stored results of zslice_preprocess, which are here
slicedir = [Tdir.pan_results,'zslice/'];

[G,S,T]=Z_get_basic_info(infile);
hh = G.h;
hh(~G.mask_rho) = -10;

var = 'salt';
switch var
    case {'salt';'temp'}
        whichgrid = 'rho';
    case {'u'}
        whichgrid = 'u';
    case {'v'}
        whichgrid = 'v';
end

s = squeeze(nc_varget(infile,'salt'));
sbot = s(1,:,:);
stop = s(end,:,:);
sfull = cat(1,sbot,s,stop);
s0 = squeeze(stop);

cvec = [30 32];
subplot(131)
pcolorcen(G.lon_rho,G.lat_rho,s0);
shading flat
caxis(cvec);
colorbar('North');
title('Surface Salinity','fontweight','bold')
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
Z_dar;
Z_info(basename,tt,T,'lr')
hold on
contour(G.lon_rho,G.lat_rho,hh,[0 0],'-k');
contour(G.lon_rho,G.lat_rho,hh,[200 200],'-k');

zlev_vec = [-50 -100];
for ii = 1:length(zlev_vec)
    zlev = zlev_vec(ii);
    load([slicedir,'zslice_',basename,'_',num2str(-zlev), ...
        '_',whichgrid,'.mat']);
    s1 = squeeze(sum(sfull.*interpmat));
    %
    subplot(1,3,ii+1)
    pcolorcen(G.lon_rho,G.lat_rho,s1);
    shading flat
    caxis(cvec);
    title([num2str(-zlev),' m Salinity'],'fontweight','bold')
    xlabel('Longitude (deg)')
    Z_dar;
    hold on
    contour(G.lon_rho,G.lat_rho,hh,[0 0],'-k');
    contour(G.lon_rho,G.lat_rho,hh,[200 200],'-k');
end






