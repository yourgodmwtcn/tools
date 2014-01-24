function plot_mask(grd,which_h,Tdir)
%plot_mask.m 4/27/2012 Parker MacCready

mask_rho = logical(grd.mask_rho);
mask_u = logical(grd.mask_u);
mask_v = logical(grd.mask_v);

switch which_h
    case 'hraw'
        h = grd.hraw;
    case 'hcarve'
        h = grd.hcarve;
    case 'h'
        h = grd.h;
    otherwise
        disp('bad which_h value')
end

h(~mask_rho) = NaN;

Z_fig;
figure

% Bathymetry
subplot(121)
Z_pcolorcen(grd.lon_rho,grd.lat_rho,h);
% note Z_pcolorcen doesn't plot tiles for the highest row and column
caxis([0 200])
colorbar
Z_dar;
Z_addcoast('combined',Tdir.coast);
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
[NR,NC] = size(grd.lon_rho);
title(['Bathymetry ',num2str(NR),' rows, ',num2str(NC),' columns']);
if 0
    hold on
    plot(grd.lon_rho(mask_rho),grd.lat_rho(mask_rho),'ok');
    plot(grd.lon_rho(~mask_rho),grd.lat_rho(~mask_rho),'xk');
    plot(grd.lon_u(mask_u),grd.lat_u(mask_u),'>k');
    plot(grd.lon_v(mask_v),grd.lat_v(mask_v),'^k');
end

% Resolution
subplot(122)
DX = 1./grd.pm; DY = 1./grd.pn;
DD = max(DX,DY);
Z_pcolorcen(grd.lon_rho,grd.lat_rho,DD);
colorbar
Z_dar;
Z_addcoast('combined',Tdir.coast);
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
title('Maximum Grid Size (m)')


