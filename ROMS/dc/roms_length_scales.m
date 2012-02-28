%   [L] = roms_length_scales(fname,varname,tindices,do_z)

function [L] = roms_length_scales(fname,varname,tindices,do_z)

if ~exist('do_z','var'), do_z = 0; warning('Ignoring z.'); end
if ~exist('tindices','var'), tindices = [1 Inf]; end

% parameters
vinfo = ncinfo(fname,varname);
dim   = length(vinfo.Size); 
slab  = 40;

% parse input
[iend,tindices,dt,nt,stride] = roms_tindices(tindices,slab,vinfo.Size(end));

time = ncread(fname,'ocean_time');
time = time([tindices(1):tindices(2)])/86400;
[xax,yax,zax,~,~] = roms_var_grid(fname,varname);

% figure out dx,dy,dz
dx = median(diff(xax));
dy = median(diff(yax));
dz = 3; % interpolate to 3m grid

L = nan(3,nt);

for i=0:iend-1
    % start and count arrays for ncread : corrected to account for stride
    [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt);
    var   = ncread(fname,varname,read_start,read_count,stride);
    
    for jj=1:size(var,4)
        ind = jj + i*slab;
        
        if do_z
            % interpolate in z
            [x,y,z] = meshgrid(xax,yax,zax);
            [xi,yi,zi] = meshgrid(xax,yax,[zax(1):dz:zax(end)]);
            txz = permute(interp3(x,y,z,permute(var(:,:,:,jj),[2 1 3]),xi,yi,zi),[2 1 3]);            
            L(3,ind) = length_scale(bsxfun(@minus,txz,mean(txz,2)),3,dz);
        else
            txz = var(:,:,:,jj);
        end
        
        L(1,ind) = length_scale(bsxfun(@minus,txz,mean(txz,2)),1,dx);
        L(2,ind) = length_scale(bsxfun(@minus,txz,mean(txz,2)),2,dy);
    end
end

% plot
figure
hold on
plot(time,L(1,:)/1000,'r');
plot(time,L(2,:)/1000,'g');
if do_z, plot(time,L(3,:),'b'); end
ylabel('Length (m/km)');
xlabel('Time (days)');
title([' Length scales : ' varname]);
legend('L_x (km)','L_y (km)','L_z (m)');