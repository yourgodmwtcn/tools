function [xax,yax,zax,xunits,yunits] = roms_var_grid(fname,varname)
    
    warning off;
    grid = roms_get_grid(fname,fname,0,1);
    
    % from John Wilkin's roms_islice.m
    % determine where on the C-grid these values lie 
    varcoords = nc_attget(fname,varname,'coordinates');
    if ~isempty(findstr(varcoords,'_u'))
      pos = 'u';
    elseif ~isempty(findstr(varcoords,'_v'))
      pos = 'v';
    elseif ~isempty(findstr(varcoords,'_rho'))
      pos = 'rho';
    else
      error('Unable to parse the coordinates variables to know where the data fall on C-grid')
    end

    switch pos
        case 'u'
            xax = grid.lon_u';
            yax = grid.lat_u';
            zax = grid.z_u;
            
            xunits = ncreadatt(fname,'x_u','units');
            yunits = ncreadatt(fname,'y_u','units'); 
            
        case 'v'
            xax = grid.lon_v';
            yax = grid.lat_v';
            zax = grid.z_v;
            
            xunits = ncreadatt(fname,'x_v','units');
            yunits = ncreadatt(fname,'y_v','units'); 

        otherwise
            xax = grid.lon_rho';
            yax = grid.lat_rho';
            
            if strcmp(varname, 'zeta');
                zax = [];
            else
                zax = grid.z_r;
            end
            xunits = ncreadatt(fname,'x_rho','units');
            yunits = ncreadatt(fname,'y_rho','units'); 
    end
    
    warning on;