% two_d_extractor.m  5//2012  Parker MacCready
%
% this is for extracting the results of ROMS simulations
% it is focused on time series of 2D fields
%

clear;
% &&&&&&&&&&& USER EDIT THIS &&&&&&&&&&&&&&
%indir = '/Users/PM/Documents/';
indir = '/pmraid1/daves/runs/'; % for skua
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

[Tdir] = pan_start; % locate things and set paths
% create output directory if needed
moordir = [Tdir.pan_results,'moorings/'];
if ~exist(moordir,'dir'); mkdir(moordir); end;

% choose run to work on and set the basename
disp('********* two_d_extractor.m *******************')
disp(' ')
if 0 % interactive version
    disp('** Select OUT where history files are **');
    pth = uigetdir(indir);
else % hardwired version
    pth = [indir,'salish_2006_4/OUT'];
end
ind = strfind(pth,'/');
basename = pth(ind(end-1)+1:ind(end)-1);
disp(' '); disp(['basename = ',basename]); disp(' ')

% get the run definition
run = roms_createRunDef('my_run',pth);
% get the year
year = datestr(run.his.nctime(1),'yyyy');
year = str2num(year);

NT = length(run.his.ncn);

% prepare output arrays
lon = run.grid.lon;
lat = run.grid.lat;
lonu = run.grid.lonu;
latu = run.grid.latu;
lonv = run.grid.lonv;
latv = run.grid.latv;
[nr,nc] = size(lon);
speed = nan(nc,nc);
max_speed = nan(nr,nc);

for tt = 4000:4720
    tt_str = ['0000',num2str(tt)];
    tt_str = tt_str(end-3:end);
    % get the correct infile
    [infile] = [run.his.dirname,run.his.basename,tt_str,'.nc'];
    if mod(tt,5)==0; disp(['  infile = ',infile]); end;
    ub = nc_varget(infile,'ubar');
    vb = nc_varget(infile,'vbar');
    ubr = interp2(lonu,latu,ub,lon,lat);
    vbr = interp2(lonv,latv,vb,lon,lat);
    speed = sqrt(ubr.^2 + vbr.*2);
    max_speed = max(max_speed,speed);
end

if 0
    figure; Z_fig;
    pcolor(lon,lat,real(max_speed));
    shading interp;
    caxis([0 1.5]);
    colorbar
    Z_dar;
    xlabel('Longitude'); ylabel('Latitude')
    title('Max depth-averaged speed (m s^{-1})');
else
    save max_speed lon lat max_speed
end
