classdef run_parent < handle
    % 12/18/2012 Parker MacCready & Sarah Giddings
    %
    % this sets the properties and methods shared by all runs
    %
    % As of 12/18/2012 I have blended in all the changes SNG made in the
    % course of creating the pnwtox runs.
    
    % initialize properties shared by all runs
    properties
        Tdir; gname; tag1; tag2; tdlims; outname; outdir;
        hmin; theta_s; theta_b; tcline; N; Vtransform; Vstretching;
        gfile; S;
    end % properties
    
    methods
        
        % note that any of the functions below can be superseded by 
        % defining a method of the same name in the child class
        
        function addInfo(rn,Tdir,gname,tag1,tag2,tdlims)
            rn.Tdir = Tdir;
            rn.gname = gname;
            rn.tag1 = tag1;
            rn.tag2 = tag2;
            rn.tdlims = tdlims;
            rn.outname = [gname,'_',tag1,'_',tag2];
        end
        
        function makeDir(rn)
            rn.outdir = [rn.Tdir.output,'mossea_run_files/', ...
                rn.outname,'/'];
            disp(' '); disp('Making run directory')
            disp(rn.outdir);
            if exist(rn.outdir)~=7; mkdir(rn.outdir); end
        end
        
        function addGrid(rn)
            grid_infile = [rn.Tdir.output,'mossea_grids/',rn.gname,'.nc'];
            rn.gfile = [rn.outdir,'grid.nc'];
            copyfile(grid_infile,rn.gfile);
        end
        
        function setScoord(rn)
            % ****** define S-coordinate parameters ***********************
            % ** THESE MUST BE IDENTICAL IN THE .IN FILE ******************
            rn.theta_s = 4;
            % THETA_S resolution goes to surface as thetaS increases
            rn.theta_b = 0.8;
            % THETA_B resolution to bed and surface equally at 1
            rn.tcline =  0;
            % TCLINE sets the depth of the pycnocline to resove
            rn.N   =     40;
            % number of vertical levels (rho)
            rn.Vtransform = 1;
            % which transformation equation to apply (1 or 2)
            rn.Vstretching = 1;
            % which stretching paramaterization to apply (1, 2, or 3). See:
            % https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
            % for transformation information.
            % Also note that if hc = 0, Vtransforms 1 and 2 are identical
            % *************************************************************
            disp(' ');
            disp(['Set vertical coordinates with N = ' num2str(rn.N) ...
                ', Ts = ' num2str(rn.theta_s) ...
                ', Tb = ' num2str(rn.theta_b)])
        end
        
        function makeS(rn)
            % make the structure "S" used by Z_s2z.m
            h = nc_varget(rn.gfile,'h'); %first get h from grid
            rn.hmin = min(h(:)); % hmin
            rn.S = Z_scoord(rn.theta_s, rn.theta_b, rn.tcline, ...
                rn.hmin, rn.N, rn.Vtransform, rn.Vstretching);
        end
        
        function makeClim(rn)
            tic
            addpath('run_functions/ocn');
            disp(' ');
            disp('Making climatology, boundary & initial files');
            indir0  = [rn.Tdir.ocn,'NCOM_global/'];
            outdir = [rn.outdir,'ocn/'];
            
            % NEED to edit to use rn.tdlims 1/7/2013
            % YUCKY
            runyear = str2num(rn.year);
            td0 = datenum(runyear,1,1,0,0,0);
            day_vec = [0 183 366];
            % will make length(day_vec)-1 clm and bry files
            subt = 1; % interval in days for subsampling
            for ii = 1:length(day_vec)-1
                day0 = day_vec(ii); day1 = day_vec(ii+1);
                tspan = [day0 day1];
                td_to_get = [day0:subt:day1] + td0;
                ts_to_get = (td_to_get-td0)*86400;
                clmname = ['ocean_clm_',num2str(ii),'.nc'];
                bryname = ['ocean_bry_',num2str(ii),'.nc'];
                disp(['  * making ',clmname]);
                make_clim(indir0, outdir, rn.gfile, ...
                    td_to_get,ts_to_get, rn.S, clmname);
                clm2bry(outdir, clmname, bryname)
                if ii==1; clmname1 = clmname; end;
            end
            % END YUCKY
            
            clm2ini(outdir,clmname1,'ocean_ini.nc')
            tclmend = toc;
            disp(['all climatology complete in ', ...
                num2str(round(tclmend)),' sec']);
        end
        
        function makeAtm(rn)
            addpath('run_functions/atm');
            disp(' ')
            disp('Making atmospheric forcing files');
            outdir = [rn.outdir,'atm/'];
            % some choices (places with year folders)
            indir0  = [rn.Tdir.atm,'mm5/tinterp_five_day/'];
            %indir0  = [rn.Tdir.atm,'mm5/tinterp_three_hour/'];
            make_atm(rn.gfile,rn.tdlims,indir0,outdir);
        end
        
        function makeTide(rn)
            addpath('run_functions/tide');
            disp(' '); disp('Making tidal forcing fields');
            indir0  = [rn.Tdir.tide,'TPXO/'];
            make_tide(rn.tdlims,rn.gfile,indir0,rn.outdir);
        end
        
        function makeRivers(rn)
            % NEED to edit to use rn.tdlims 1/7/2013
            addpath('run_functions/river');
            disp(' '); disp('Adding rivers');
            riv_dir  = [rn.Tdir.river];
            riverFile = [riv_dir,'riverShapes.mat'];
            flowFile = [riv_dir,'riverFlow_1998_2009.mat'];
            disp('use file from grid')
            tspan = [datenum(str2double(rn.year),1,1)
                datenum(str2double(rn.year)+1,1,1)];
            %make rivers without initializing/ramping flow
            make_rivers(tspan,rn.year,rn.S,rn.gfile,rn.outdir, ...
                riverFile,flowFile,0)
            %alters IC inside river areas
            riverpolyFile = [riv_dir,'riverPolygons.mat'];
            ininame = 'ocn/ocean_ini.nc';
            ini_rivers(rn.year,rn.S,rn.outdir,ininame,flowFile, ...
                riverFile,riverpolyFile)
        end
        
        function addDye(rn)
            % do nothing
            disp(' '); disp('no dye added')
        end        
        
    end % methods
    
end % classdef