function [Fb_v] = Z_make_Fb_v(infile)
% 12/16/2011 Parker MacCready
%
% this makes a lat,lon array of turbulent buoyancy flux per unit volume
% (units W m-3) for a given infile

[G,S,T] = Z_get_basic_info(infile);
zbot = -G.h;
zeta = nc_varget(infile,'zeta');
z_rho = roms_z(G.h,zeta,S.Cs_r);

salt =  nc_varget(infile,'salt');
temp =  nc_varget(infile,'temp');
% get the potential density ("temp" is potential temperature in ROMS)
[rho] = Z_make_potdens(salt,temp);

% form drhodz
DZ = diff(z_rho); % on the w grid (no top and bottom)
Drho = diff(rho); % same
drhodz = Drho ./ DZ;
% then get the diffusivity
K = nc_varget(infile,'AKs');
% form Fb
g = 9.8;
Fb_v = -g * K(2:end-1,:,:) .* drhodz; % should be Watts/m3

