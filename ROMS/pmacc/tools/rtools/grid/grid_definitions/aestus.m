classdef aestus < grid_parent & grid_utilities
    % 5/21/2012 Parker MacCready
    %
    % This is an idealized estuary
    
    methods
        % just paste in any of the functions from the method section
        % of grid_parent.m, edit as needed, and the version here will take
        % precedence
        
        
        function addRiverFile(grd)
            grd.riverFile = 'none_used';
        end
        
        function addTopoFiles(grd)
            grd.topoFiles = 'none_used';
        end
        
        function addMaskFile(grd)
            grd.maskFile = [];
        end
        
        function makeLonLat(grd)
            grd.hmin = 4;
            grd.spherical = 'T';
            lon1 = -1; lon2 = 2; lat1 = 44; lat2 = 46;
            delta = 3000; % grid spacing (m)
            RE = Z_RE(lat1,lat2); % earth radius (m)
            dlat = delta / (pi*RE/180);
            dlon = dlat / cos(pi*0.5*(lat1+lat2)/180);
            lonr = [lon1:dlon:lon2];
            latr = [lat1:dlat:lat2];
            [grd.lon_rho,grd.lat_rho] = meshgrid(lonr,latr);
        end
        
        function makeHraw(grd)
            hraw = grd.hmin*ones(size(grd.lon_rho));
            % carve the estuary channel
            hraw(grd.lat_rho>=44.9 & grd.lat_rho<=45.1) = 20;
            % carve the continental slope
            hraw(grd.lon_rho<0) = -200*grd.lon_rho(grd.lon_rho<0);
            % save results
            grd.hraw = hraw;
            grd.hcarve = grd.hraw; % in case we don't carve rivers
            grd.h = grd.hraw; % in case we don't smooth
            grd.h(grd.h<grd.hmin) = grd.hmin; % do this before smoothing
            grd.h(isnan(grd.h)) = grd.hmin; %fill in nan's on
            % eastern edge of grid where there is no ocean (needed?)
        end
        
        function carveRivers(grd)
            disp('NOT carving rivers')
        end
        
        function editMask(grd)
            disp('NOT editing the mask')
        end

                
    end % methods
    
end % classdef