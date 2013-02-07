%       [xax,yax,zax,tax,xunits,yunits] = roms_var_grid(fname,varname)
% returns grids as matrices of appropriate size for variable

function [xax,yax,zax,tax,xunits,yunits] = dc_roms_var_grid(fname,varname,tindex)
    
    warning off;
    
    if ~exist('tindex','var'), tindex = 0; end
    
    if strcmp(varname,'pv')
        pos = 'p';
    elseif strcmp(varname,'eke') || strcmp(varname,'mke') || strcmp(varname,'pe')
        pos = 'e';
    elseif strcmp(varname,'vor')
        pos = 'q';
    else
        grid = roms_get_grid(fname,fname,tindex,1);

        % from John Wilkin's roms_islice.m
        % determine where on the C-grid these values lie 
        varcoords = nc_attget(fname,varname,'coordinates');
        if ~isempty(strfind(varcoords,'_u'))
          pos = 'u';
        elseif ~isempty(strfind(varcoords,'_v'))
          pos = 'v';
        elseif ~isempty(strfind(varcoords,'_w'))
          pos = 'w';
        elseif ~isempty(strfind(varcoords,'_rho'))
          pos = 'r'; % rho
        else
          error('Unable to parse the coordinates variables to know where the data fall on C-grid')
        end
    end
    
    N = size(grid.z_r,1);

    switch pos
        case 'u'
            xax = repmat(grid.lon_u',[1 1 N]);
            yax = repmat(grid.lat_u',[1 1 N]);
            zax = permute(grid.z_u,[3 2 1]);
            tax = ncread(fname,'ocean_time');
            
            try
                xunits = ncreadatt(fname,'lon_u','units');
                yunits = ncreadatt(fname,'lat_u','units'); 
            catch ME
                xunits = ncreadatt(fname,'x_u','units');
                yunits = ncreadatt(fname,'y_u','units'); 
            end
            
        case 'v'
            xax = repmat(grid.lon_v',[1 1 N]);
            yax = repmat(grid.lat_v',[1 1 N]);
            zax = permute(grid.z_v,[3 2 1]);
            tax = ncread(fname,'ocean_time');
            
            try
                xunits = ncreadatt(fname,'lon_v','units');
                yunits = ncreadatt(fname,'lat_v','units'); 
            catch ME
                xunits = ncreadatt(fname,'x_v','units');
                yunits = ncreadatt(fname,'y_v','units');
            end                
            
        case 'w'
            xax = repmat(grid.lon_rho',[1 1 N+1]);
            yax = repmat(grid.lat_rho',[1 1 N+1]);
            zax = permute(grid.z_w,[3 2 1]);
            tax = ncread(fname,'ocean_time');
            
            try
                xunits = ncreadatt(fname,'lon_rho','units');
                yunits = ncreadatt(fname,'lat_rho','units'); 
            catch ME
                xunits = ncreadatt(fname,'x_rho','units');
                yunits = ncreadatt(fname,'y_rho','units'); 
            end

        case 'r'
            
            tax = ncread(fname,'ocean_time');
            
            if strcmp(varname, 'zeta');
                xax = grid.lon_rho';
                yax = grid.lat_rho';
                zax = [];
            else
                xax = repmat(grid.lon_rho',[1 1 N]);
                yax = repmat(grid.lat_rho',[1 1 N]);
                zax = permute(grid.z_r,[3 2 1]);
            end
            
            try
                xunits = ncreadatt(fname,'lon_rho','units');
                yunits = ncreadatt(fname,'lat_rho','units'); 
            catch ME
                xunits = ncreadatt(fname,'x_rho','units');
                yunits = ncreadatt(fname,'y_rho','units'); 
            end
        
        case 'p'
            xax = ncread(fname,'x_pv');
            yax = ncread(fname,'y_pv');
            zax = ncread(fname,'z_pv');
            tax = ncread(fname,'ocean_time');
            
            xunits = ncreadatt(fname,'x_pv','units');
            yunits = ncreadatt(fname,'y_pv','units');   
            
       case 'e'
            xax = ncread(fname,'x_en');
            yax = ncread(fname,'y_en');
            zax = ncread(fname,'z_en');
            tax = ncread(fname,'t_en');
            xunits = ncreadatt(fname,'x_en','units');
            yunits = ncreadatt(fname,'y_en','units');   
            
       case 'q'
            xax = ncread(fname,'xvor');
            yax = ncread(fname,'yvor');
            zax = ncread(fname,'zvor');
            tax = ncread(fname,'ocean_time');
            xunits = ncreadatt(fname,'xvor','units');
            yunits = ncreadatt(fname,'yvor','units');  
    end
    
    xunits = strrep(xunits,'_',' ');    
    yunits = strrep(yunits,'_',' ');  
    
    warning on;