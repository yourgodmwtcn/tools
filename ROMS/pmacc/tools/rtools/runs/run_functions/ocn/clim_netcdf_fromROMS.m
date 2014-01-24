function [V,AAlon,AAlat,AAmask,AAlon_old,AAlat_old,AAmask_old,k] = ...
    clim_netcdf_fromROMS(varname,newgrid,G)
% 12/6/2011  Parker MacCready
%
% returns variable-specific information
% SNG 19 June 2012 updated to work for extracting from a ROMS file

% Read in the new grid matrices
lon_rho = nc_varget(newgrid,'lon_rho');
lat_rho = nc_varget(newgrid,'lat_rho');
mask_rho = nc_varget(newgrid,'mask_rho');

lon_u = nc_varget(newgrid,'lon_u');
lat_u = nc_varget(newgrid,'lat_u');
mask_u = nc_varget(newgrid,'mask_u');

lon_v = nc_varget(newgrid,'lon_v');
lat_v = nc_varget(newgrid,'lat_v');
mask_v = nc_varget(newgrid,'mask_v');

% Read in the old grid matrices
lon_rho_old = G.lon_rho;
lat_rho_old = G.lat_rho;
mask_rho_old = logical(G.mask_rho);

lon_u_old = G.lon_u;
lat_u_old = G.lat_u;
mask_u_old = logical(G.mask_u);

lon_v_old = G.lon_v;
lat_v_old = G.lat_v;
mask_v_old = logical(G.mask_v);

do_depth = 1; % flat to indicate that there is depth
do_addtimevar = 1; % include the time variable (set to 0 for v and vbar)

switch varname
    case 'salt'
        invarname = 'salt';
        ncvarname = 'salt';
        nclongname = 'salinity climatology';
        ncunits = 'PSU';
        nctimename = 'salt_time';
        ncxiname = 'xi_rho'; ncetaname = 'eta_rho';
        AAlon = lon_rho; AAlat = lat_rho; AAmask = mask_rho;
        AAlon_old = lon_rho_old; AAlat_old = lat_rho_old; AAmask_old = mask_rho_old;
    case 'temp'
        invarname = 'temp';
        ncvarname = 'temp';
        nclongname = 'potential temperature climatology';
        ncunits = 'Celsius';
        nctimename = 'temp_time';
        ncxiname = 'xi_rho'; ncetaname = 'eta_rho';
        AAlon = lon_rho; AAlat = lat_rho; AAmask = mask_rho;
        AAlon_old = lon_rho_old; AAlat_old = lat_rho_old; AAmask_old = mask_rho_old;
    case 'u'
        invarname = 'u';
        ncvarname = 'u';
        nclongname = 'u-momentum component climatology';
        ncunits = 'meter second-1';
        nctimename = 'v3d_time';
        ncxiname = 'xi_u'; ncetaname = 'eta_u';
        AAlon = lon_u; AAlat = lat_u; AAmask = mask_u;
        AAlon_old = lon_u_old; AAlat_old = lat_u_old; AAmask_old = mask_u_old;
    case 'v'
        invarname = 'v';
        ncvarname = 'v';
        nclongname = 'v-momentum component climatology';
        ncunits = 'meter second-1';
        nctimename = 'v3d_time';
        do_addtimevar = 0;
        ncxiname = 'xi_v'; ncetaname = 'eta_v';
        AAlon = lon_v; AAlat = lat_v; AAmask = mask_v;
        AAlon_old = lon_v_old; AAlat_old = lat_v_old; AAmask_old = mask_v_old;
    case 'zeta'
        invarname = 'zeta';
        ncvarname = 'zeta';
        nclongname = 'sea surface height climatology';
        ncunits = 'meter';
        nctimename = 'zeta_time';
        do_depth = 0; % 2D variable
        ncxiname = 'xi_rho'; ncetaname = 'eta_rho';
        AAlon = lon_rho; AAlat = lat_rho; AAmask = mask_rho;
        AAlon_old = lon_rho_old; AAlat_old = lat_rho_old; AAmask_old = mask_rho_old;
    case 'ubar'
        invarname = 'ubar';
        ncvarname = 'ubar';
        nclongname = 'vertically averaged u-momentum climatology';
        ncunits = 'meter second-1';
        nctimename = 'v2d_time';
        do_depth = 0; % 2D variable
        ncxiname = 'xi_u'; ncetaname = 'eta_u';
        AAlon = lon_u; AAlat = lat_u; AAmask = mask_u;
        AAlon_old = lon_u_old; AAlat_old = lat_u_old; AAmask_old = mask_u_old;
    case 'vbar'
        invarname = 'vbar';
        ncvarname = 'vbar';
        nclongname = 'vertically averaged v-momentum climatology';
        ncunits = 'meter second-1';
        nctimename = 'v2d_time';
        do_depth = 0; % 2D variable
        do_addtimevar = 0;
        ncxiname = 'xi_v'; ncetaname = 'eta_v';
        AAlon = lon_v; AAlat = lat_v; AAmask = mask_v;
        AAlon_old = lon_v_old; AAlat_old = lat_v_old; AAmask_old = mask_v_old;
end

switch varname
    case  {'salt', 'temp', 'u', 'v', 'zeta'}
        % fill out k which will be used for nearest neighbor extrap
        k = dsearchn([AAlon_old(AAmask_old) AAlat_old(AAmask_old)], ...
            [AAlon_old(~AAmask_old) AAlat_old(~AAmask_old)]);
    otherwise
        k = [];
end

% pack the netcdf info in a structure "V"
V.invarname = invarname;
V.ncvarname = ncvarname; V.nclongname = nclongname; 
V.ncunits = ncunits; V.nctimename = nctimename;
V.do_depth = do_depth; V.do_addtimevar = do_addtimevar;
V.ncxiname = ncxiname; V.ncetaname = ncetaname;
