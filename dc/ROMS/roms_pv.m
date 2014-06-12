% calculates Ertel PV at interior RHO points (horizontal plane) and midway between rho points in the vertical
%       [pv] = roms_pv(fname,tindices)

function [pv,xpv,ypv,zpv] = roms_pv(fname, tindices, outname, ftype)

if ~exist('ftype', 'var'), ftype = 'avg'; end
if isdir(fname)
    dirname = fname;
    fnames = roms_find_file(dirname, ftype);
    fname = [dirname '/' fnames{1}];
    tpv = dc_roms_read_data(dirname,'ocean_time', [], {}, [], [], ftype);
    dirflag = 1;
else
    % extract dirname
    index = strfind(fname,'/');
    dirname = fname(1:index(end));
    tpv = dc_roms_read_data(fname,'ocean_time');
    dirflag = 0;
end
vinfo = ncinfo(fname,'u');
s     = vinfo.Size;
dim   = length(s); 
slab  = roms_slab(fname,0) + 6;

warning off
grid = roms_get_grid(fname,fname,1,1);
warning on

% parse input
if ~exist('tindices','var'), tindices = []; end

[iend,tindices,dt,~,stride] = roms_tindices(tindices,slab,length(tpv));

rho0  = ncread(fname,'R0');
tpv = tpv(tindices(1):tindices(2));
f   = ncread(fname,'f',[1 1],[Inf Inf]);

xname = 'x_pv'; yname = 'y_pv'; zname = 'z_pv'; tname = 'ocean_time';
xrname = 'x_rv'; yrname = 'y_rv'; zrname = 'z_rv';

xdname = 'xpv'; ydname = 'ypv'; zdname = 'zpv'; tdname = 'tpv';
xdrname = 'xrv'; ydrname = 'yrv'; zdrname = 'zrv';

grid1.xv = repmat(grid.x_v',[1 1 grid.N]);
grid1.yv = repmat(grid.y_v',[1 1 grid.N]);
grid1.zv = permute(grid.z_v,[3 2 1]);

grid1.xu = repmat(grid.x_u',[1 1 grid.N]);
grid1.yu = repmat(grid.y_u',[1 1 grid.N]);
grid1.zu = permute(grid.z_u,[3 2 1]);

grid1.xr = repmat(grid.x_rho',[1 1 grid.N]);
grid1.yr = repmat(grid.y_rho',[1 1 grid.N]);
grid1.zr = permute(grid.z_r,[3 2 1]);

grid1.zw = grid.z_w;
grid1.s_w = grid.s_w;
grid1.s_rho = grid.s_rho;

totvol = sum(grid.dV(:));
%% setup netcdf file

if ~exist('outname','var') || isempty(outname), outname = 'ocean_vor.nc'; end
outname = [dirname '/' outname];

if exist(outname,'file')
    %in = input('File exists. Do you want to overwrite (1/0)? ');
    in = 1;
    if in == 1, delete(outname); end
end

% better packing = better compression?
%pv_scale = 1e-11;
%rv_scale = 1e-4;
% steal chunking from ROMS
rhoinfo = ncinfo(fname, 'rho');
chunksize = rhoinfo.ChunkSize;

nccreate(outname,'pv', 'Format','netcdf4', 'DeflateLevel',5,'Shuffle',true,...
    'Dimensions', {xdname s(1)-1 ydname s(2)-2 zdname s(3)-1 tdname length(tpv)}, ...
    'ChunkSize', chunksize, 'Datatype', 'single');
nccreate(outname,'rv', 'Format','netcdf4', 'DeflateLevel',5,'Shuffle',true,...
    'Dimensions', {xdrname s(1) ydrname s(2)-1 zdrname s(3)-1 tdname length(tpv)}, ...
    'ChunkSize', chunksize, 'Datatype', 'single');
nccreate(outname,xname,'Dimensions',{xdname s(1)-1 ydname s(2)-2 zdname s(3)-1});
nccreate(outname,yname,'Dimensions',{xdname s(1)-1 ydname s(2)-2 zdname s(3)-1});
nccreate(outname,zname,'Dimensions',{xdname s(1)-1 ydname s(2)-2 zdname s(3)-1});
nccreate(outname,xrname,'Dimensions',{xdrname s(1)  ydrname s(2)-1 zdrname s(3)-1});
nccreate(outname,yrname,'Dimensions',{xdrname s(1)  ydrname s(2)-1 zdrname s(3)-1});
nccreate(outname,zrname,'Dimensions',{xdrname s(1)  ydrname s(2)-1 zdrname s(3)-1});
nccreate(outname,tname,'Dimensions',{tdname length(tpv)});
nccreate(outname,'intPV','Dimensions',{tdname length(tpv)});

ncwriteatt(outname,'pv','Description','Ertel PV calculated from ROMS output');
ncwriteatt(outname,'pv','coordinates',[xname ' ' yname ' ' zname ' ' tname]);
ncwriteatt(outname,'pv','units','N/A');
%ncwriteatt(outname,'pv','scale_factor',pv_scale);
ncwriteatt(outname,'rv','Description','Relative voritcity, vx-uy');
ncwriteatt(outname,'pv','coordinates',['x_rv y_rv z_rv ocean_time']);
ncwriteatt(outname,'rv','units','1/s');
%ncwriteatt(outname,'rv','scale_factor',rv_scale);
ncwriteatt(outname,xname,'units',ncreadatt(fname,'x_u','units'));
ncwriteatt(outname,yname,'units',ncreadatt(fname,'y_u','units'));
ncwriteatt(outname,zname,'units','m');
ncwriteatt(outname,tname,'units','s');
ncwriteatt(outname,'intPV','Description', ...
    'time series of volume averaged PV over entire domain.');
fprintf('\n Created file : %s\n', outname);

%% calculate pv

misc = roms_load_misc(fname);
ticstart = tic;
for i=0:iend-1
    disp(['i = ' num2str(i) '/' num2str(iend-1)]);
    [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt);
    tstart = read_start(end);
    tend   = read_start(end) + read_count(end) -1;
    
    if dirflag
        u = dc_roms_read_data(dirname,'u',[tstart tend],{},[],grid, ...
                              ftype, 'single');
        v = dc_roms_read_data(dirname,'v',[tstart tend],{},[],grid, ...
                              ftype, 'single');
        
        try
            rho = dc_roms_read_data(dirname,'rho',[tstart tend],{},[],grid, ...
                                    ftype, 'single');
        catch ME
            rho = rho0 -rho0 * misc.Tcoef* ...
                dc_roms_read_data(dirname,'temp',[tstart tend],{},[],grid, ...
                                  ftype, 'single');
        end
    else
        u = ncread(fname,'u',read_start,read_count,stride);
        v = ncread(fname,'v',read_start,read_count,stride);
        try
            rho = ncread(fname,'rho',read_start,read_count,stride); % theta
        catch ME
            rho = -misc.Tcoef*ncread(fname,'temp',read_start,read_count,stride);
        end
    end

    disp('Calculating pv...');
    pvstart = tic;
    [pv,xpv,ypv,zpv,rvor] = pv_cgrid(grid1,u,v,rho,f,rho0);
    toc(pvstart);
    
    if i == 0
        % write grid
        ncwrite(outname,xname,xpv);
        ncwrite(outname,yname,ypv);
        ncwrite(outname,zname,zpv);
        ncwrite(outname,'ocean_time',tpv);
    end

    disp('Writing data...');
    ticstartwrt = tic;
    ncwrite(outname,'pv',pv,read_start);
    ncwrite(outname,'rv',rvor,read_start);
    toc(ticstartwrt);
    
    pvdV = bsxfun(@times,pv,avg1(grid.dV(2:end-1,2:end-1,:),3));
    
    intPV(tstart:tend) = squeeze(nansum(nansum(nansum(pvdV,1),2),3))./totvol; 
end
ncwrite(outname,'intPV',intPV);
toc(ticstart);
%save pv.mat pv xpv ypv zpv tpv intPV
fprintf('\n Wrote file : %s \n\n',outname);

    %% old code
    
%     pv1    = avgx(avgz(bsxfun(@plus,avgy(vx - uy),f)))  .*  (tz(2:end-1,2:end-1,:,:));
%     pv2    = (-1)*;
%     pv3    = uz.*avgz(tx);
    %pv = double((pv1 + avgy(pv2(2:end-1,:,:,:)) + avgx(pv3(:,2:end-1,:,:)))./avgz(lambda(2:end-1,2:end-1,:,:))); 