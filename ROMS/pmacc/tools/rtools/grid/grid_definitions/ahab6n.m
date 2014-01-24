classdef ahab6n < grid_parent & grid_utilities
    % 4/27/2012 Neil Banas & Parker MacCready
    
    methods
                
        function addTopoFiles(grd)
            grd.topoFiles = {'psdem/TTP_Regional_27m.mat', ...
                'psdem/PS_91m.mat', ...
                'cascadia/cascadia_gridded.mat', ...
                'smith_sandwell/pnw_smithsand.mat'};
        end
        
        function addMaskFile(grd)
            grd.maskFile = 'ahab6n_goodmask.nc';
        end

        function makeLonLat(grd)
            
            res = 0.6;
            
            lonmin = -124;
			lonmax = -122.1;
			latmin = 46.9;
			latmax = 49.5;
			x =       [lonmin  -123.15  -122.3  lonmax];
			resx_km = [     2      res     res       1];
			y =       [latmin   47.1   48.1   48.7  latmax];
			resy_km = [     2    res    res      1       2];
			r = 1.1;
            
            km2m = 1000; %convert km to m
            RE = Z_RE(latmin,latmax); %earth radius at mean lat
            latm = (latmin+latmax)/2; %mean latitude
            %dl = 1/80;
            grd.hmin = 4;
            grd.spherical = 'T';
            
            % make the lon-lat grid
            resx = resx_km ./ (pi*RE/km2m/180) ./ cosd(latm);
            resy = resy_km ./ (pi*RE/km2m/180);
            lonr = resolution2grid(x,resx,r);
            latr = resolution2grid(y,resy,r);
            
            [grd.lon_rho,grd.lat_rho] = meshgrid(lonr,latr);

			% estimate speed relative to ptx_high
			demand = prod(size(grd.lon_rho)) / (570*249) / ...
                min( min(resy_km), min(resx_km) );
			disp(['ahab6n: run time = ' num2str(demand) ...
                ' relative to ptx_high']);
			% estimate speed relative to ptx_highT
			demand = prod(size(grd.lon_rho)) / (381*174) / ...
                min( min(resy_km), min(resx_km) );
			disp(['ahab6n: run time = ' num2str(demand) ...
                ' relative to ptx_highT']);
        end
                
    end % methods
    
end % classdef
