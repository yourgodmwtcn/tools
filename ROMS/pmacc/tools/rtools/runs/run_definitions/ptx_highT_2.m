classdef ptx_highT_2 < run_parent
    % 12/18/2012 Parker MacCready
    %
    % this s mostly obsolete, but contains the code to add dye, as an
    % example
    
    methods
         
        function setScoord(rn)
            rn.theta_s = 4;
            rn.theta_b = 0.8;
            rn.tcline =  0;
            rn.N   =     27;
            rn.Vtransform = 1;
            rn.Vstretching = 1;
        end  
        
        function addDye(rn)
            tic
            addpath('run_functions/dye/');
            %DyeRiverNames lists the names of the rivers you want to add dye to
            DyeRiverNames = {'columbia','fraser','others'}; %add dye to columbia, fraser rivers independently and the rest together
            %NESW is a vector indicating dye on (1) or off (0) the north, east, south, west boundaries respectively.
            NESW = [0 0 1 1]; %dye on south and west boundaries
            %get string list of all rivers included
            riv_dir  = [rn.Tdir.river];
            riverFile = [riv_dir,'ps_2005_riverFile_ptxhigh3.mat'];
            load(riverFile)
            for ri = 1:length(rivers)
                AllRiverNames{ri} = rivers(ri).name;
            end
            clear rivers
            numdye = length(DyeRiverNames) + nansum(NESW); %total dyes
            ocean_dir = [rn.outdir,'ocn/'];
            %find how many clm files are in ocean_dir
            d = dir(ocean_dir);
            for i = 1:length(d)
                if length(d(i).name)>=9
                    ni(i) = strcmp(d(i).name(1:9),'ocean_clm');
                else
                    ni(i) = 0;
                end
            end
            numclm = nansum(ni);
            %loop through each clm and bry file
            for ci = 1:numclm
                %clm file will have zero dye everywhere
                clmfile = ([ocean_dir 'ocean_clm_' num2str(ci) '.nc']);
                add_dyes_to_clm(clmfile,NESW,numdye)
                %bry file will have ones on the desired boundaries only
                bryfile = ([ocean_dir 'ocean_bry_' num2str(ci) '.nc']);
                add_dyes_to_bry(bryfile,NESW,numdye);
            end
            %the ini file will have zero dye everywhere
            inifile = [ocean_dir,'ocean_ini.nc'];
            add_dyes_to_ini(inifile,5);
            %the river file will have dye in the designated rivers but not the brys
            river_file = [rn.outdir,'rivers.nc'];
            add_dyes_to_rivers(river_file,AllRiverNames,DyeRiverNames,numdye)
            tdyeend = toc;
            disp(['dye addition complete in ',num2str(round(tdyeend)),' sec']);
        end        

                
    end % methods
    
end % classdef