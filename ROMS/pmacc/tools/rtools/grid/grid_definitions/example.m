classdef example < grid_parent & grid_utilities
    % 4/27/2012 Parker MacCready
    %
    % This is a template to use when creating your own grid specification
    % file.
    
    methods        
        % just paste in any of the functions from the method section
        % of grid_parent.m, edit as needed, and the version here will take
        % precedence
        
        function makeLonLat(grd)
            grd.hmin = 4;
            grd.spherical = 'T';
            lon1 = -127; lon2 = -122; lat1 = 45; lat2 = 50;
            delta = 5000; % grid spacing (m)
            RE = Z_RE(lat1,lat2); % earth radius (m)
            dlat = delta / (pi*RE/180);
            dlon = dlat / cos(pi*0.5*(lat1+lat2)/180);
            lonr = [lon1:dlon:lon2];
            latr = [lat1:dlat:lat2];
            [grd.lon_rho,grd.lat_rho] = meshgrid(lonr,latr);
        end
        
    end % methods
    
end % classdef