classdef example_1 < run_parent
    % 4/27/2012 Parker MacCready
    
    methods
        
        function setScoord(rn)
            rn.theta_s = 4;
            rn.theta_b = 0.8;
            rn.tcline =  0;
            rn.N   =     27;
            rn.Vtransform = 1;
            rn.Vstretching = 1;
        end
        
        function makeClim(rn)
            addpath('run_functions/ocn');
            disp('Making climatology, boundary & initial files');
            indir  = [rn.Tdir.ocn,'NCOM_global/pro',rn.year,'/'];
            runyear = str2num(rn.year);
            day_vec = [0:20:40]; % debugging
            subt = 10; % interval in days for subsampling
            for ii = 1:length(day_vec)-1
                day0 = day_vec(ii); day1 = day_vec(ii+1);
                tspan = [day0 day1];
                clmname = ['ocean_clm_',num2str(ii),'.nc'];
                bryname = ['ocean_bry_',num2str(ii),'.nc'];
                disp(['  * making ',clmname]);
                make_clim(indir, rn.outdir, rn.gfile, ...
                    tspan, subt, rn.S, clmname);
                clm2bry(rn.outdir, clmname, bryname)
                if ii==1; clmname1 = clmname; end;
            end
            clm2ini(rn.outdir,clmname1,'ocean_ini.nc')
        end
        
        
    end % methods
    
end % classdef