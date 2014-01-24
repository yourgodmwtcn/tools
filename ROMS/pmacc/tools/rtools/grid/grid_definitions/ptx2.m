classdef ptx2 < grid_parent & grid_utilities
    % 10/28/2011 SNG ptxhigh Turbo grid
    % (1.5x coarser resolution than ptxhigh)
    % Vertical resolution is also 1.5x coarser and thus rx0 can be 1.5x
    % higher (i.e. rx0 = 0.15)
    % some edits 4/27/2012 PM for consistency with grid_maker.m
    
    methods
        
        function addRiverFile(grd)
            grd.riverFile = 'ps_2005_riverFile_ptxhigh3.mat';
        end
        
        function addMaskFile(grd)
            grd.maskFile = 'ptx_highT_goodmask.nc';
        end

        function makeLonLat(grd)
            disp('Making a stretched grid')
            % higher res ps grid, try swath approach instead of parabola
            %set grid limits
            lonmin = -127.4;
            lonmax = -122;
            latmin = 43;
            latmax = 50;
            km2m = 1000; %convert km to m
            RE = Z_RE(latmin,latmax); %earth radius at mean lat
            latm = (latmin+latmax)/2; %mean latitude
            grd.hmin = 4;
            grd.spherical = 'T';
            %set min and max resolutions (km) in each tanh function
            %they are going to be exactly 1.5x the ptxhigh values
            minresw = 1.5;
            maxresw = 4.5;
            minrese = minresw;
            maxrese = 3;
            minress = 1.5;
            maxress = 3;
            minresn = minress;
            maxresn = 3;
            %to set up the tanh profiles, need to calculate delta_res
            dresw = maxresw-minresw;
            drese = maxrese-minrese;
            dress = maxress-minress;
            dresn = maxresn-minresn;
            %set scale factors to determine the width over which the
            %tanh decays (higher number = wider tanh)
            %widens it to about 2*scalew km on either side, for example
            %scalew = 40 makes a tanh nearly 200km wide, scalew = 20
            %makes a tanh nearly 100km wide
            scalew = 40;
            scalee = 30;
            scales = 50;
            scalen = 40;
            % x is km with zero at lower left corner
            % y is km with zero at lower left corner
            lonr_0 = linspace(lonmin,lonmax,200);
            latr_0 = linspace(latmin,latmax,200);
            x=(pi*RE*cosd(latm)/180)*(lonr_0-(lonmin))/km2m;
            y=(pi*RE/180)*(latr_0-(latmin))/km2m;
            xt1 = x(dsearchn(lonr_0',-125.7));
            %target lon to start tanh decay at
            yt1 = y(dsearchn(latr_0',45.5));
            xt2 = x(dsearchn(lonr_0',-123.25));
            %target lon to end tanh decay at
            yt2 = y(dsearchn(latr_0',49));
            % shape is a function to multiply dl by
            shape_w = ((minresw+dresw/2)+(dresw/2).*tanh(-(x-xt1)/scalew));
            shape_e = ((minrese+drese/2)+(drese/2).*tanh((x-xt2)/scalee));
            shape = max(shape_e,shape_w);
            %
            i = 1;
            lonr(i) = lonmin;
            while lonr(i)<=lonmax
                i = i+1;
                shape_fact = interp1(lonr_0,shape,lonr(i-1));
                lonr(i)=lonr(i-1)+km2m*shape_fact*(180/(pi*RE*cosd(latm)));
            end
            
            % Get the latitude for an isotropic grid
            shape_s = ((minress+dress/2)+(dress/2).*tanh(-(y-yt1)/scales));
            shape_n = ((minresn+dresn/2)+(dresn/2).*tanh((y-yt2)/scalen));
            shape = max(shape_n,shape_s);
            
            i=1;
            latr(i)=latmin;
            while latr(i)<=latmax
                i=i+1;
                shape_fact = interp1(latr_0,shape,latr(i-1));
                latr(i) = latr(i-1)+km2m*shape_fact*(180/(pi*RE));
            end
            
            [grd.lon_rho,grd.lat_rho] = meshgrid(lonr,latr);
            
        end
        
    end % methods
    
end % classdef