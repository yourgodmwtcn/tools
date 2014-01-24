classdef grid_parent < handle
    % 10/27/2011 Parker MacCready
    %
    % this sets the properties and methods shared by all grids
    
    % initialize properties shared by all grids
    properties
        gname; Tdir; hraw; hcarve; h; spherical; f;
        pm; pn; hmin; xl; el;
        lon_rho; lat_rho; mask_rho;
        lon_u; lat_u; mask_u; lon_v; lat_v; mask_v;
        lon_psi; lat_psi; mask_psi;
        x_rho; y_rho; x_u; y_u; x_v; y_v; x_psi; y_psi;
        carve_done = 'F'; smooth_done = 'F'; mask_edited = 'F';
        riverFile; topoFiles; maskFile;
    end % properties
    
    methods
        
        % note that any of the functions below can be superseded by
        % defining methods of the same name in the child class
        
        function addRiverFile(grd)
            grd.riverFile = 'riverShapes.mat';
            % this has info on river channels to carve
        end
        
        function addTopoFiles(grd)
            grd.topoFiles = {'cascadia/cascadia_gridded.mat',...
                'smith_sandwell/pnw_smithsand.mat'};
            % can be a cell array or single string
            % make sure to list these in order of preference, best first
        end
        
        function addMaskFile(grd)
            grd.maskFile = [];
            % copy this method to your child and name a grid file (.nc)
            % to reuse that mask
        end

        function makeLonLat(grd)
            disp('Making grid')
            grd.hmin = 4;
            grd.spherical = 'T';
            lon1 = -127; lon2 = -122; lat1 = 43; lat2 = 50;
            delta = 5000; % grid spacing (m)
            RE = Z_RE(lat1,lat2); % earth radius (m)
            dlat = delta / (pi*RE/180);
            dlon = dlat / cos(pi*0.5*(lat1+lat2)/180);
            lonr = [lon1:dlon:lon2];
            latr = [lat1:dlat:lat2];
            [grd.lon_rho,grd.lat_rho] = meshgrid(lonr,latr);
        end
        
        function hraw = makeHraw_oneSource(grd, matfile)
            % interpolates bathymetry from one gridded source file.
            % doesn't save! that's what makeHraw() is for.
            % specify source bathymetry either as an absolute file reference
            % or relative to the 'topo' directory specified in toolstart.m
            if exist(matfile)
                topofile = matfile;
            else
                topofile = [grd.Tdir.topo matfile];
            end
            topo = load(topofile);
            hraw = interp2(topo.lon,topo.lat,-topo.z, ...
                grd.lon_rho,grd.lat_rho);
        end
        
        function makeHraw(grd)
            % makes raw bathymetry from one or more gridded source files.
            if ~iscell(grd.topoFiles)
                % single bathymetry source
                hraw = grd.makeHraw_oneSource(grd.topoFiles);
            else
                % start with bathymetry source #1
                hraw = grd.makeHraw_oneSource(grd.topoFiles{1});
                for i=2:length(grd.topoFiles)
                    % fill in area that remains uncovered
                    % with the next bathymetry source
                    hraw2 = grd.makeHraw_oneSource(grd.topoFiles{i});
                    hraw(isnan(hraw)) = hraw2(isnan(hraw));
                end
            end
            % save results
            grd.hraw = hraw;
            grd.hcarve = grd.hraw; % in case we don't carve rivers
            grd.h = grd.hraw; % in case we don't smooth
            grd.h(grd.h<grd.hmin) = grd.hmin; % do this before smoothing
            grd.h(isnan(grd.h)) = grd.hmin; %fill in nan's on
            % eastern edge of grid where there is no ocean (needed?)
        end
        
        function makeMask(grd)
            disp('Making intial mask')
            hh = grd.hraw; % depth (m)
            hmin = grd.hmin;
            grd.mask_rho = double(hh > hmin); % mask_rho = 1 in water
            [grd.mask_u,grd.mask_v,grd.mask_psi] = ...
                Z_mask_uvp(grd.mask_rho);
        end
        
        function carveRivers(grd)
            disp('carving rivers')
            river_file = [grd.Tdir.river grd.riverFile];
            %use grd.h which at this point is equivalent to hraw but all
            %values are >= hmin
            [rout, mask_rho, hcarve] = Z_carve_river_channels(grd.h, ...
                grd.lon_rho,grd.lat_rho,grd.mask_rho,river_file);
            grd.hcarve = hcarve;
            grd.mask_rho = double(mask_rho);
            [grd.mask_u,grd.mask_v,grd.mask_psi] = ...
                Z_mask_uvp(grd.mask_rho);
            grd.carve_done = 'T';
        end
        
        function editMask(grd)
            grd.writeNetCDF; % editmask_rtools needs this
            disp('edit the mask')
            grid_file = [grd.Tdir.output 'mossea_grids/' grd.gname,'.nc'];
            coast_file = [grd.Tdir.coast 'pnw_coast_combined.mat'];
            disp('Press RETURN after exiting editmask')
            editmask_rtools(grid_file, coast_file);
            pause
            grd.mask_edited = 'T';
            grd.mask_rho = nc_varget(grid_file,'mask_rho');
            [grd.mask_u,grd.mask_v,grd.mask_psi] = ...
                Z_mask_uvp(grd.mask_rho);
        end
        
        function loadMask(grd)
            ncname = [grd.Tdir.output 'mossea_grids/' grd.maskFile];
            disp(['loading mask from ' ncname]);
            grd.mask_edited = 'T';
            grd.mask_rho = nc_varget(ncname,'mask_rho');
            [grd.mask_u,grd.mask_v,grd.mask_psi] = ...
                Z_mask_uvp(grd.mask_rho);
        end
        
        function smooth(grd)
            %set the max r factor to smooth to
            rfact = 0.15;
            disp(['Smoothing the grid using r = ' num2str(rfact)])
            %compute areas at each rho point
            areaMatrix = (1./grd.pm).*(1./grd.pn);
            %GRID_PlusMinusScheme_rx0 is part of the LP_Bathymetry package
            grd.h = GRID_PlusMinusScheme_rx0(grd.mask_rho, ...
                grd.hcarve, rfact, areaMatrix);
            grd.smooth_done = 'T';
        end
        
    end % methods
    
end % classdef