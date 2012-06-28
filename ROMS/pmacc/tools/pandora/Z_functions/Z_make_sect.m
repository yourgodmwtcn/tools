function [sect] = Z_make_sect(infile,lon,lat,var);
% 4/27/2010  Parker MacCready
%
% makes a section of a variable along a specified path which may be
% irregular; uses SNCTOOLS
%
% INPUT:
% infile [text] the name of the netcdf input file (ROMS history file)
% lon [vector double]
% lat [vector double of the same length]
% var [text] the name of a 3D variable
%
% OUTPUT is a structure:
% sect.time_datenum [scalar] the datenum format time
% sect.varname [string] = name of the variable
% sect.infile [string] = name of the input file
% sect.lon [matrix double] lon vector, now at all depths
% sect.lat [matrix double] lat vector, now at all depths
% sect.dist [matrix double] along-section distance matrix (m)
%   measured along the lat,lon path provided as input
% sect.ang_rad [matrix double] of the angle of the section measured as
%   radians relative to East (positive counter clockwise)
% sect.z [matrix double] z postions of all points (m, positive up)
%   including the free surface and bottom
%
% sect.var [matrix double] the variable on the section, with the
%   variable name being one of a subset of ROMS history variable names
%   If the variable name is 'uv' then the elements of the array are complex
%   with u being the real part and v being the imaginary part
% thus the along-section and normal-to-section velocities would be calculated as
%   ca = cos(sect.ang_rad); sa = sin(sect.ang_rad);
%   us = ca.*real(sect.var) + sa.*imag(sect.var);
%   un = ca.*imag(sect.var) - sa.*real(sect.var);
%
% 5/26/2010 now using the faster "roms_z" function

% get basic grid info
[G,S,T] = Z_get_basic_info(infile);

sect.time_datenum = T.time_datenum;
sect.varname = var;
sect.infile = infile;

% bathymetry, SSH, and z levels
h = nc_varget(infile,'h');
zbot = -h;
zeta = nc_varget(infile,'zeta');
if 1
    z_rho = roms_z(h,zeta,S.Cs_r);
    z_w = roms_z(h,zeta,S.Cs_w);
else
    [z_rho,z_w] = Z_s2z_mat(h,zeta,S);
end
zzbot(1,:,:) = zbot;
zzeta(1,:,:) = zeta;
z_full = cat(1,zzbot,z_rho,zzeta);

% distance and angle vectors
clat = cos(2*pi*mean(lat)/360);
RE = 6371e3; % radius of Earth (m)
dlon = diff(lon); dlat = diff(lat);
dx = clat*RE*2*pi*dlon/360; dy = RE*2*pi*dlat/360;
delta = sqrt(dx.^2 + dy.^2);
dist = [0, cumsum(delta)];
ang_rad = atan2(dy,dx);
% extrapolate for the angle of the end point
ang_rad = [ang_rad ang_rad(end)];

% interpolate on the section for different classes of variables
switch var
    case {'salt','temp','rho'}
        c = nc_varget(infile,var);
        for ii = 1:S.N
            this_c = squeeze(c(ii,:,:));
            this_c(~G.mask_rho) = NaN;
            sect_c(ii,:) = interp2(G.lon_rho,G.lat_rho,this_c,lon,lat);
        end
        sect.var = cat(1,sect_c(1,:),sect_c,sect_c(end,:));
        sect.lon = repmat(lon,S.N+2,1);
        sect.lat = repmat(lat,S.N+2,1);
        sect.dist = repmat(dist,S.N+2,1);
        sect.ang_rad = repmat(ang_rad,S.N+2,1);
        for ii = 1:S.N+2
            this_z = squeeze(z_full(ii,:,:));
            sect.z(ii,:) = interp2(G.lon_rho,G.lat_rho,this_z,lon,lat);
        end
    case {'uv'}
        c = nc_varget(infile,'u');
        for ii = 1:S.N
            this_c = squeeze(c(ii,:,:));
            this_c(~G.mask_u) = NaN;
            sect_u(ii,:) = interp2(G.lon_u,G.lat_u,this_c,lon,lat);
        end
        c = nc_varget(infile,'v');
        for ii = 1:S.N
            this_c = squeeze(c(ii,:,:));
            this_c(~G.mask_v) = NaN;
            sect_v(ii,:) = interp2(G.lon_v,G.lat_v,this_c,lon,lat);
        end
        sect.var = sect_u + i*sect_v;
        % note that we set zero velocity at the bottom and no shear at the
        % top (neither of which are strictly correct, but they are fine for
        % plotting)
        sect.var = cat(1,0*sect.var(1,:),sect.var,sect.var(end,:));
        sect.lon = repmat(lon,S.N+2,1);
        sect.lat = repmat(lat,S.N+2,1);
        sect.dist = repmat(dist,S.N+2,1);
        sect.ang_rad = repmat(ang_rad,S.N+2,1);
        for ii = 1:S.N+2
            this_z = squeeze(z_full(ii,:,:));
            sect.z(ii,:) = interp2(G.lon_rho,G.lat_rho,this_z,lon,lat);
        end
    case {'w','AKs','AKv'}
        c = nc_varget(infile,var);
        for ii = 1:S.N+1
            this_c = squeeze(c(ii,:,:));
            this_c(~G.mask_rho) = NaN;
            sect.var(ii,:) = interp2(G.lon_rho,G.lat_rho,this_c,lon,lat);
        end
        sect.lon = repmat(lon,S.N+1,1);
        sect.lat = repmat(lat,S.N+1,1);
        sect.dist = repmat(dist,S.N+1,1);
        sect.ang_rad = repmat(ang_rad,S.N+1,1);
        for ii = 1:S.N+1
            this_z = squeeze(z_w(ii,:,:));
            sect.z(ii,:) = interp2(G.lon_rho,G.lat_rho,this_z,lon,lat);
        end
    case {'Fb'} % the buoyancy flux - a derived variable (will be on z_w)
        % first get the density
        salt =  nc_varget(infile,'salt');
        temp =  nc_varget(infile,'temp');
        [rho] = Z_make_potdens(salt,temp); % use potential density
        % note that "temp" actually is theta already 2/11/2012
        for ii = 1:S.N
            this_rho = squeeze(rho(ii,:,:));
            this_rho(~G.mask_rho) = NaN;
            sect_rho(ii,:) = interp2(G.lon_rho,G.lat_rho,this_rho,lon,lat);
        end
        % get z on rho
        for ii = 1:S.N
            this_z = squeeze(z_rho(ii,:,:));
            sect_zr(ii,:) = interp2(G.lon_rho,G.lat_rho,this_z,lon,lat);
        end
        % form drhodz
        DZ = diff(sect_zr); % on the w grid (no top and bottom)
        Drho = diff(sect_rho); % same
        sect_drhodz = Drho ./ DZ;
        sect_drhodz_full = cat(1,0*sect_drhodz(1,:),sect_drhodz,sect_drhodz(end,:));
        % then get the diffusivity
        K = nc_varget(infile,'AKs');
        for ii = 1:S.N+1
            this_K = squeeze(K(ii,:,:));
            this_K(~G.mask_rho) = NaN;
            sect_K(ii,:) = interp2(G.lon_rho,G.lat_rho,this_K,lon,lat);
        end
        % form Fb
        g = 9.8;
        sect.var = -g * sect_K .* sect_drhodz_full; % should be Watts/m3
        % etc
        sect.lon = repmat(lon,S.N+1,1);
        sect.lat = repmat(lat,S.N+1,1);
        sect.dist = repmat(dist,S.N+1,1);
        sect.ang_rad = repmat(ang_rad,S.N+1,1);
        for ii = 1:S.N+1
            this_z = squeeze(z_w(ii,:,:));
            sect.z(ii,:) = interp2(G.lon_rho,G.lat_rho,this_z,lon,lat);
        end
end


