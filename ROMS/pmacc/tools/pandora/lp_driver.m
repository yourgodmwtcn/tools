% lp_driver.m 12/14/2011 Parker MacCready
%
% makes a tidally-averaged version of ROMS output files

clear;
% &&&&&&&&&&& USER EDIT THIS &&&&&&&&&&&&&&
indir = '/Users/PM/Documents/Salish/runs/';
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

[Tdir] = pan_start; % locate things and set paths
Tdir.pan_lp = [Tdir.pan_top,'lp/'];
if ~exist(Tdir.pan_lp,'dir'); mkdir(Tdir.pan_lp); end;

% choose run to work on and set the basename
disp('********* lp_driver.m *******************'); disp(' ')
if 1 % interactive version
    disp('** Select OUT where history files are **');
    pth = uigetdir(indir);
else % hardwired version
    pth = [indir,'salish_2006_3/OUT'];
end
ind = strfind(pth,'/');
basename = pth(ind(end-1)+1:ind(end)-1);
disp(' '); disp(['basename = ',basename]); disp(' ')

nickname = basename;
run = roms_createRunDef(nickname,pth);

outdir = [Tdir.pan_lp,basename,'/OUT/'];
if exist(outdir)==7; rmdir(outdir,'s'); end;
mkdir(outdir)

%filterWindow = 3/24; % for testing
filterWindow = 40/24; % a standard 40 hour tidal filter

suffix = 'lp';
run.hislp = roms_hanning(run.his,filterWindow,[],suffix,outdir);

