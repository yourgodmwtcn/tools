% mooring_extractor_batch.m
%
%this code does mooring extractios using Z_get_moor.m for
%all runs in indirlist and for all mooring locations in mlocations and
%saves them as matfiles in your outpath
%sng 6 June 2012, similar to ROMS_mooring_extration_driver, but now calls
%Z_get_moor.m, PMs improved mooring extractor. 
%
%note that the user needs to change the input directories (indirlist)
%and the output directory(outdir)
%the user also needs to choose the mooring locations either from a
%predefined list in mooringLocations.m or their own lat/lon locations (see
%the examples below)
%Users also need to ensure that their toolstart.m includes the right paths
%
%to run this on skua in the background go to its directory and run:
%matlab<ROMS_mooring_extraction_driver_Zgetmoor.m>&mextract_log&
%
% renamed by PM 6/8/2012 to fit into pandora

%% User input
%set list of input directories to loop through
%note the list can contain more than one run, for example:
% indirlist = {'/pmraid4/sarahgid/runs/ptx_highT_2_2005/OUT', ...
%      '/pmraid3/sarahgid/runs/ptx_highT_2_2005_addxxobc/OUT';};
indirlist = {'/pmraid4/sarahgid/runs/ptx_highT_2_2005/OUT', ...
    '/pmraid3/sarahgid/runs/ptx_highT_2_2005_addxxobc/OUT';};

%set list of mooring locations to extract
%moorings already included in mooringLocations.m 
%(EH1, EH2, EH3, EH4, RN, RC, RS, J2C, A1, racerocks)
%alterntively a user can enter their own lat/lon locations and mooring names
%example pre-defined mooring location input
% mlocations = {'EH1','EH2','EH3','EH4','RN','RC','RS','J2C','A1','racerocks'};
% OR a shorter list: mlocations = {'EH2','EH4'};
%example user input locations (i.e. mooring locations not defined in
%mooringLocations.m) Note that alternatively, mooringLocations.m can be updated
% mlocations = {'S300_100','S100_100','S300_300','S100_300','W47_100','W47_300'};
% lat = [43.1 43.1 43.27 43.27 47 47];
% lon = [-124.8582 -124.5643 -124.71 -124.5579 -127.27 -127];
% OR another example: mlocations = {'scidac1','nemo_new',ai0'};
% lon = [-125.9 -125 -122.7];
% lat = [48.5 48 48.13]; 
mlocations = {'EH1','EH2','EH3','EH4','RN','RC','RS','J2C','A1','racerocks'};
    
%set output directory
outpath = '/pmraid1/sarahgid/runs/PostProcess/ROMS_moorings/';

%% set the appropriate paths
% set global paths and directories
addpath('../alpha/'); 
Tdir = toolstart;
addpath(Tdir.tools);

%% do the mooring extractions
%loop through input directories
for inum = 1:length(indirlist)
    indir = indirlist{inum};
    
    %loop through locations to extract
    for mnum = 1:length(mlocations)
        mloc = mlocations{mnum};
        
        %get run information
        % assumes that "indir" is something like:
        % /Users/PM/Documents/Salish/runs/ptx_med_2005_1/OUT
        % then the lines below find the string right before "OUT", which is the
        % basename, the main identifier of the run
        ind = strfind(indir,'/');
        basename = indir(ind(end-1)+1:ind(end)-1);
        disp(' '); disp(['basename = ',basename]);
        disp(' '); disp(['mooring = ',mloc]); disp(' ')
        % also figure out if we are extracting low-passed files
        oname = indir(ind(end)+1:end);
        if strcmp(oname,'OUT_lp') || strcmp(oname,'OUT_lp_hanning')
            % NOTE the assumed naming convention
            time_flag = 'ave';
            hisbasename = 'ocean_his_lp_';
        else
            time_flag = 'his';
            hisbasename = 'ocean_his_';
        end
        % get the run definition
        runi = roms_createRunDef('my_run',indir,hisbasename);
        % get the year of the run
        yeari = datestr(runi.his.nctime(1),'yyyy');
        yeari = str2double(yeari);

        %get the mooring locations
        run Z_mooringLocations;
        
        %run ROMS_mooring_extractor
        [M] = Z_get_moor(runi,lon_moor,lat_moor);
        % pack the results in a structure
        M.run = runi;
        M.basename = basename;
        % rename for backward consistency
        mod_moor = M;
        
        %save the output mooring file
        if strcmp(time_flag,'ave');
            save([outpath,basename,'_lp_',mloc,'.mat'],'mod_moor');
        else
            save([outpath,basename,'_',mloc,'.mat'],'mod_moor');
        end
    end
end
