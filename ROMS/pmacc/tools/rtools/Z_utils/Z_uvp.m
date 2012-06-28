function [xu,yu,xv,yv,xp,yp]=Z_uvp(xr,yr)

% 10/17/2008  Parker MacCready
% transforms coordinates from the rho-grid to the u-, v-, and psi-grids

xu = xr(:,1:end-1) + diff(xr,1,2)/2;
yu = yr(:,1:end-1);
xv = xr(1:end-1,:);
yv = yr(1:end-1,:) + diff(yr,1,1)/2;
xp = xu(1:end-1,:);
yp = yv(:,1:end-1);
