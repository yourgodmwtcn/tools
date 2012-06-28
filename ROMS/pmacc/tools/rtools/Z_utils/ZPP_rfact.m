function [rx0, rx1] = ZPP_rfact(h, mask_rho , S)
%-------------------------------
%  [rx0,rx1]=rfact(h, S)
%
% h = grid bathymetry, mask_rho = mask grid
% S = grid parameter structure from info.mat
% ***note, if you just want rx0, omit S use rx0 = ZPP_rfact(h, mask_rho)
%
% fxn to calculate the grid "stiffnesses" rx0 and rx1 from the bathymetry
% and vertical structure fxn used
%
% rx0 is like deltaH / H and ROMS suggests to be between 0-0.4
% rx1 is related to vertical structure (Haney #) and ROMS suggests 3-7
%
%-------------------------------

do_rx1 = 1;
if nargin < 3  %if don't have S, can still calculate rx0
    do_rx1 = 0;    
end

%first get sizes of grid
[M,L] = size(h);
Mm  = M-1;
Mmm = M-2;
Lm  = L-1;
Lmm = L-2;
[mask_u,mask_v,mask_psi] = Z_mask_uvp(mask_rho);

% then calculate r factors, use h, in stiffness.F they use bottom z_w value
%---- my_rx0=MAX(my_rx0,ABS((z_w(i,j,0)-z_w(i-1,j,0))/            &
%     &                            (z_w(i,j,0)+z_w(i-1,j,0))))  --------
rx = abs(h(1:M,2:L)-h(1:M,1:Lm))./(h(1:M,2:L)+h(1:M,1:Lm));
rx = rx.*mask_u;
ry = abs(h(2:M,1:L)-h(1:Mm,1:L))./(h(2:M,1:L)+h(1:Mm,1:L));
ry = ry.*mask_v;
%
rx0 = max(max(rx(:)),max(ry(:)));

%---------------------------------------------------------
if(do_rx1)
%next, calculate rx1
[z_rho,z_w] = Z_s2z_mat(h,0.*h,S);
rxx1 = 0; ryy1 = 0; 
for k = 2:S.N+1
    zz0 = squeeze(z_w(k-1,:,:));
    zz1 = squeeze(z_w(k,:,:));
    rxtemp = abs((zz1(1:M,2:L) - zz1(1:M,1:Lm) + zz0(1:M,2:L) - zz0(1:M,1:Lm)) ./ ...
                                 (zz1(1:M,2:L) + zz1(1:M,1:Lm) - zz0(1:M,2:L) - zz0(1:M,1:Lm)));
    rxx1 = max(rxx1, max(max(rxtemp.*mask_u)));
    rytemp = abs((zz1(2:M,1:L) - zz1(1:Mm,1:L) + zz0(2:M,1:L) - zz0(1:Mm,1:L)) ./ ...
                                 (zz1(2:M,1:L) + zz1(1:Mm,1:L) - zz0(2:M,1:L) - zz0(1:Mm,1:L)));
    ryy1 = max(ryy1, max(max(rytemp.*mask_v)));
end
%
rx1 = max(rxx1, ryy1);
end
% from Utility/stiffness.F
% ------- my_rxx1=MAX(my_rx1,ABS((z_w(i,j,k  )-z_w(i-1,j,k  )+       &
%         &                           z_w(i,j,k-1)-z_w(i-1,j,k-1))/      &
%         &                           (z_w(i,j,k  )+z_w(i-1,j,k  )-       &
%         &                           z_w(i,j,k-1)-z_w(i-1,j,k-1)))) ------