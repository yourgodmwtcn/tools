classdef nesttest_1 < run_parent
    % 21 June 2012 SNG
    
    methods
         
        function setScoord(rn)
            % ****** define S-coordinate parameters **************************
            % ** THESE MUST BE IDENTICAL IN THE .IN FILE *********************
            rn.theta_s = 4;
            % THETA_S resolution goes to surface as thetaS increases
            rn.theta_b = 0.8;
            % THETA_B resolution to bed and surface equally at 1
            rn.tcline =  0;
            % TCLINE sets the depth of the pycnocline to resove
            rn.N   =     20;
            % number of vertical levels (rho)
            rn.Vtransform = 1;
            % which transformation equation to apply (1 or 2)
            rn.Vstretching = 1;
            % which stretching paramaterization to apply (1, 2, or 3). See:
            % https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
            % for transformation information.
            % Also note that if hc = 0, Vtransforms 1 and 2 are identical
            % ****************************************************************
            
            % make the structure "S" used by Z_s2z.m
            h = nc_varget(rn.gfile,'h'); %first get h from grid
            rn.hmin = min(h(:)); % hmin
            rn.S = Z_scoord(rn.theta_s, rn.theta_b, rn.tcline, ...
                rn.hmin, rn.N, rn.Vtransform, rn.Vstretching);
        end
        
        function makeClim(rn)
            %this is different from the run_parent.m makeClim because it
            %uses ROMS output files as the boundary and initial conditions
            %rather than a global model (NCOM)
            addpath('run_functions/ocn');
            disp('Making climatology, boundary & initial files');
            %put the directory path for your roms run here
            roms_dir = '/pmraid3/sarahgid/runs/ptx_highT_2_2005_addxxobc/OUT';
            %desired output directory, do not change this!
            outdir = [rn.outdir,'ocn/'];
            %input desired time span here (matlab datenum)
            tspan = [datenum('jan 1 2005') datenum('dec 31 2005')];
            %input subsampling time (hrs)
            subt = 720;
            %climatology, boundary files
            clmname = 'ocean_clm.nc';
            bryname = 'ocean_bry.nc';
            disp(['  * making ',clmname]);
            %use make_clim_fromROMS_new to export data from the ROMS model
            %to be used at the boundaries and initial conditions and
            %climatology files
            make_clim_fromROMS_new(roms_dir, outdir, clmname, ...
                rn.gfile, rn.S, tspan, subt);
            %make boundary file from climatology file
            clm2bry(outdir, clmname, bryname)
            %make initial condition file from the climatolgy file
            clm2ini(outdir,clmname,'ocean_ini.nc')
        end
        
        function makeAtm(rn)
            %input the atmospheric forcing files
            addpath('run_functions/atm');
            outdir = [rn.outdir,'atm/'];
            disp('Making atmospheric forcing files');
            %choose 5 day or 3 hr time steps
            dth  = 120;
            switch dth
                case 120
                    indir00 = 'tinterp_five_day/';
                case 3
                    indir00 = 'tinterp_three_hour/';
            end
            indir0  = [rn.Tdir.atm,'mm5/',indir00]; % full of year folders
            %make atmospheric forcing files
            make_atm(rn.gfile,rn.year,dth,indir0,outdir);
        end
        
        function makeTide(rn)
            disp('not making tides because using ssh from roms output');
        end
        
        function makeRivers(rn)
            %only do this if you have rivers in your domain!
            addpath('run_functions/river');
            disp('Adding rivers');
            riv_dir  = [rn.Tdir.river];
            %this river file must just contain the rivers you plan to
            %include
            riverFile = [riv_dir,'ps_2005_riverFile_nesttest.mat'];
            disp('use file from grid')
            make_rivers(rn.year,rn.S,rn.gfile,riv_dir,rn.outdir,riverFile,0);
            %alters IC inside river areas
            riverpolyFile = [riv_dir,'PSriverPolys_ptxhigh3.mat'];
            ininame = 'ocn/ocean_ini.nc';
            Z_PS_ini_rivers(rn.year, rn.S, rn.outdir, ininame, riv_dir, riverFile, riverpolyFile)
        end
        
        function addDye(rn)  
            disp('not adding any dye')
        end
        
    end % methods
    
end % classdef