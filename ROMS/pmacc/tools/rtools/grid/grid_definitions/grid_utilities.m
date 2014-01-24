classdef grid_utilities < handle
    % 4/27/2012 Parker MacCready
    %
    % methods that should never have to be changed
    
    methods
        
        function addGname(grd,gname)
            grd.gname = gname;
        end
        
        function addTdir(grd,Tdir)
            grd.Tdir = Tdir;
        end
        
        function makeGridExtras(grd)
            if(grd.spherical);
                lon_rho = grd.lon_rho; lat_rho = grd.lat_rho;
                lat1 = min(lat_rho(:,1)); lat2 = max(lat_rho(:,1));
                lon1 = min(lon_rho(1,:)); lon2 = max(lon_rho(1,:));
                RE = Z_RE(lat1,lat2);
                cl=cos(pi*0.5*(lat1+lat2)/180);
                grd.xl = pi*cl*RE*(lon2-lon1)/180; % domain x length (m)
                grd.el = pi*RE*(lat2-lat1)/180; % domain y length (m)
                x_rho = pi*RE*cl*(lon_rho-lon_rho(1,1))/180;
                y_rho = pi*RE*(lat_rho-lat_rho(1,1))/180;
                grd.x_rho = x_rho; grd.y_rho = y_rho;
                % make uvp coordinates
                [grd.lon_u,grd.lat_u,grd.lon_v,grd.lat_v, ...
                    grd.lon_psi,grd.lat_psi] = ...
                    Z_uvp(lon_rho,lat_rho);
                [grd.x_u,grd.y_u,grd.x_v,grd.y_v, ...
                    grd.x_psi,grd.y_psi] = Z_uvp(x_rho,y_rho);
                % make metrics
                [grd.pm,grd.pn] = Z_pmpn(grd.x_u,grd.y_v);
                grd.f = sw_f(lat_rho);
            else
                disp('Non-spherical not implemented yet')
                % reorder for analytical grid (assume we created x_rho and
                % y_rho and then want to construct lon and lat).
            end
        end
        
        function writeNetCDF(grd)
            ncfile_out = [grd.Tdir.output, ...
                'mossea_grids/',grd.gname,'.nc'];
            disp(['Creating ' ncfile_out]);
            my_mode = bitor ( nc_clobber_mode, nc_64bit_offset_mode );
            nc_create_empty ( ncfile_out, my_mode );
            % add dimensions
            nx = length(grd.lon_rho(1,:));
            ny = length(grd.lat_rho(:,1));
            nc_add_dimension( ncfile_out, 'xi_rho', nx);
            nc_add_dimension( ncfile_out, 'eta_rho', ny);
            nc_add_dimension( ncfile_out, 'xi_u', nx-1);
            nc_add_dimension( ncfile_out, 'eta_u', ny);
            nc_add_dimension( ncfile_out, 'xi_v', nx);
            nc_add_dimension( ncfile_out, 'eta_v', ny-1);
            nc_add_dimension( ncfile_out, 'xi_psi', nx-1);
            nc_add_dimension( ncfile_out, 'eta_psi', ny-1);
            % add the variables and write to them
            varlist = {'spherical','xl','el', ...
                'hraw', 'hcarve', 'h', ...
                'f','pm','pn', ...
                'x_rho','y_rho','x_u','y_u','x_v','y_v','x_psi','y_psi', ...
                'lon_rho','lat_rho','lon_u','lat_u', ...
                'lon_v','lat_v','lon_psi','lat_psi', ...
                'mask_rho','mask_u','mask_v','mask_psi'};
            for ii = 1:length(varlist)
                varname = varlist{ii};
                varstruct = Z_make_varstruct(varname);
                nc_addvar(ncfile_out,varstruct);
                nc_varput(ncfile_out, varname, eval(['grd.',varname]));
            end
        end
        
    end % methods
    
end % classdef