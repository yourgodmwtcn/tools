function [L] = roms_length_scales(fname,varname,tindices)

    
% parameters
vinfo = ncinfo(fname,varname);
dim   = length(vinfo.Size); 
slab  = 40;

%warning off
%grid = roms_get_grid(fname,fname,0,1);
%warning on

% parse input
[iend,tindices,dt,nt,stride] = roms_tindices(tindices,slab,vinfo.Size(end));

time = ncread(fname,'ocean_time');
time = time([tindices(1):tindices(2)])/86400;
[xax,yax,zax,~,~] = roms_var_grid(fname,varname);

% figure out dx,dy,dz
dx = median(diff(xax));
dy = median(diff(yax));
dz = 3; % interpolate to 3m grid

for i=0:iend-1
    % start and count arrays for ncread : corrected to account for stride
    [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt);
    var   = ncread(fname,varname,read_start,read_count,stride);
    
    % interpolate in z
    
    for jj=1:size(var,4)
        L(1,jj) = length_scale(var(:,:,:,i),1,dx);
        L(2,jj) = length_scale(var(:,:,:,i),2,dy);
        L(3,jj) = length_scale(var(:,:,:,i),3,dz);
    end
end