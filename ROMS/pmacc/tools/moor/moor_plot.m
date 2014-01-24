% moor_plot.m 2/5/2013 Parker MacCready
%
% plots the results of a mooring extraction

clear;
%close all;
moor_start;

% select the mooring file
[fn,pth]=uigetfile([Tdir.moor_out,'*.mat'],'Select mooring file...');
infile = [pth,fn];
load(infile);
% do this for backward compatibility
if exist('mod_moor','var'); M = mod_moor; clear mod_moor; end;

% select the plot file
[fn_p,pth_p]=uigetfile(['Z_plot_code/*.m'],'Select plot code...');
plot_file = strrep(fn_p,'.m',''); % used in an "eval" call below

% make the plot(s)
nplots = length(M);
for mmm = 1:nplots
    MM = M(mmm);
    MM.run = M(1).run;
    MM.basename = M(1).basename;
    eval([plot_file,'(MM);']);
end

% code below needs to be reconfigured for multiple moorings
%
% save the plot if desired
if 0 % don't do this unless you want to
    if nplots == 1
        plotit = input('Save a plot? [1 = save, RETURN = do not save]');
        if isempty(plotit); plotit = 0; end;
        if plotit
            set(gcf,'PaperPositionMode','auto');
            print('-djpeg100',[Tdir.moor_out,fn,'_',plot_file,'.jpg']);
        end
    end
end


