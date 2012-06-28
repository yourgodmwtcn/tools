function [mask_u,mask_v,mask_psi] = Z_mask_uvp(mask_rho);

% 10/17/2008  Parker MacCready
% This makes the u_, v_, and psi_mask from rho_mask

mask_u = mask_rho(:,1:end-1).*mask_rho(:,2:end);
mask_v = mask_rho(1:end-1,:).*mask_rho(2:end,:);
mask_psi = mask_u(1:end-1,:).*mask_u(2:end,:);
