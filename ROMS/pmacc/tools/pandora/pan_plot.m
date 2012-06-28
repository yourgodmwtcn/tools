% pan_plot.m  12/7/2011  Parker MacCready
%
% a general tool for selecting a ROMS history file (or group of files) and
% plotting it using a veriety of different plotting codes
%
% if more than one history file is chosen then a sqeuence of plots is saved
% in a new directory, and these may be combined into a movie using other
% software
%
% it assumes that there is one history file per time step

clear
%close all

% set paths, get Tdir structure of directories, and
% create output directories if needed
% (relies on tools/aplha/toolstart.m)
Tdir = pan_start;

% choose file or files to plot
disp('********* pan_plot.m ********************************')
disp(' ')
indir = 'E:\Work\CattlePass\runs\';

% choose which run to use, and set the basename
[fn,pth]=uigetfile([indir,'*.nc'], ...
    'Select NetCDF file or files...','multiselect','on');
% ASSUMES that "pth" is something like:
% /Users/PM/Documents/Salish/runs/ptx_med_2005_1/OUT/
% then the lines below find the string right before "OUT", which is the
% basename, the main identifier of the run
ind = strfind(pth,'\');
basename = pth(ind(end-2)+1:ind(end-1)-1);
disp(' '); disp(['basename = ',basename]); disp(' ')
% also figure out if we are plotting low-passed files
if isempty(strfind(fn,'lp')) % NOTE the assumed naming convention
    time_flag = '';
else
    time_flag = '_lp';
end

% if you selected more than one file it will "make a movie" which means
% creating a directory full of separate jpeg's which one can then plot
% using something like "open image sequence" in QuickTime 7
make_movie = 0;
if iscell(fn); make_movie = 1; end

% choose which plotting code to use
[fn_p,pth_p] = uigetfile('Salish_plot_code/*.m','Select Plotting code...');
addpath(pth_p);
plot_file = strrep(fn_p,'.m',''); % used in an "eval" call below

% default initialization of plot properties
Z_fig(10); figure;
set(gcf,'position',[100 100 1000 600]);

% determine how many files to plot
if make_movie; ntt = size(fn,2); else; ntt = 1; end;

for ii = 1:ntt % MOVIE loop start (or just make single plot)
    if make_movie;
        infile = [pth,fn{ii}];
    else
        infile = [pth,fn];
    end
    
    tt = str2num(infile(end-6:end-3));
    
    % make the plot
    switch plot_file
        % *** EDITED TO NEW STANDARDS *********** 
        case 'basic_mooring'
            % requires that you have created mooring files, and that you
            % have river files
            if ii==1;
                % get the river file
                T = Z_get_time_info(infile);
                year = datestr(T.time_datenum,'yyyy'); % a character string
                riv = load([Tdir.river,'PS_flow_', ...
                    year,'.mat']);
                % get a mooring file
                [fn_m,pth_m] = uigetfile([Tdir.pan_results,'moorings/', ...
                    basename,'*.mat'],'Select mooring file...');
                load([pth_m,fn_m]); % gets mod_moor structure
                % reset figure size
                Z_fig(12);
                set(gcf,'position',[50 50 1300 700]);
            end;
            eval([plot_file,'(Tdir,infile,basename,tt,riv,mod_moor);']);
        otherwise
            % this is the default plotting call, using the plot_file that
            % you selected
            eval([plot_file,'(Tdir,infile,basename,tt);']);
    end
    
    % PRINTING
    %
    % NOTE: printing with -dtiff instead of -djpeg makes a cleaner plot,
    % but the file is about 4x as big.  Maybe best for publication.
    %
    % NOTE: there is a weird bug in that the printed figure is increased in
    % size by a factor 5/3 over what is on the sceen (as determined by
    % pixel count).
    %
    if make_movie == 0
        % choose to save a copy of the plot
        plotit = input('Save a plot? [1 = save, RETURN = do not save]');
        if isempty(plotit); plotit = 0; end;
        if plotit
            set(gcf,'PaperPositionMode','auto');
            print('-dtiff',[Tdir.pan_fig,plot_file,'_',basename, ...
                time_flag,'_',num2str(tt),'.tif']);
        end
    else % make a folder of jpegs for a movie
        if ii==1
            outdir = [Tdir.pan_mov,plot_file,'_',basename,time_flag];
            if exist(outdir)==7; rmdir(outdir,'s'); end;
            mkdir(outdir);
        end
        set(gcf,'PaperPositionMode','auto');
        print('-djpeg90',[outdir,'/',num2str(tt),'.jpg']);
        if ii<length(fn); clf; end;
    end
    
end % MOVIE loop end
disp('DONE');

