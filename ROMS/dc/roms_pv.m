% calculates Ertel PV
%       [pv] = roms_pv(fname,tindices)

function [pv] = roms_pv(fname,tindices)

% parameters
lam = 'rho';
vinfo = ncinfo(fname,'u');
s     = vinfo.Size;
dim   = length(s); 
slab  = 20;

warning off
grid = roms_get_grid(fname,fname,0,1);
warning on

% parse input
if ~exist('tindices','var'), tindices = []; end

[iend,tindices,dt,nt,stride] = roms_tindices(tindices,slab,vinfo.Size(end));

%%
% caps indicates domain integrated values
EKE = nan(nt,1);
MKE = EKE;
PE  = EKE;

R0   = ncread(fname,'R0');
time = ncread(fname,'ocean_time');
time = time([tindices(1):tindices(2)])/86400;
f    = ncread(fname,'f',[1 1],[Inf Inf]);
f = mean(f(:));

pv = nan([s(1)-1 s(2)-2 s(3)-1 tindices(2)-tindices(1)+1]);

for i=0:iend-1
    [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt);
    tstart = read_start(end);
    tend   = read_start(end) + read_count(end) -1;
    
    u      = ncread(fname,'u',read_start,read_count,stride);
    v      = ncread(fname,'v',read_start,read_count,stride);
    lambda = ncread(fname,lam,read_start,read_count,stride); % theta

    % calculate gradients
    vx    = bsxfun(@rdivide,diff(v,1,1),diff(grid.x_v',1,1)); %diff(v,1,1)./repmat(diff(grid.x_v',1,1),[1 1 s(3) s(4)]);
    vy    = bsxfun(@rdivide,diff(v,1,2),diff(grid.y_v',1,2)); %diff(v,1,2)./repmat(diff(grid.y_v',1,2),[1 1 s(3) s(4)]);
    vz    = bsxfun(@rdivide,diff(v,1,3),permute(diff(grid.z_v,1,1),[3 2 1])); %diff(v,1,3)./repmat(permute(diff(grid.z_v,1,1),[3 2 1]),[1 1 1 s(4)]);

    ux    = bsxfun(@rdivide,diff(u,1,1),diff(grid.x_u',1,1)); %diff(v,1,1)./repmat(diff(grid.x_v',1,1),[1 1 s(3) s(4)]);
    uy    = bsxfun(@rdivide,diff(u,1,2),diff(grid.y_u',1,2)); %diff(v,1,2)./repmat(diff(grid.y_v',1,2),[1 1 s(3) s(4)]);
    uz    = bsxfun(@rdivide,diff(u,1,3),permute(diff(grid.z_u,1,1),[3 2 1])); %diff(v,1,3)./repmat(permute(diff(grid.z_v,1,1),[3 2 1]),[1 1 1 s(4)]);

    tx    = bsxfun(@rdivide,diff(lambda,1,1),diff(grid.x_rho',1,1)); %diff(v,1,1)./repmat(diff(grid.x_v',1,1),[1 1 s(3) s(4)]);
    ty    = bsxfun(@rdivide,diff(lambda,1,2),diff(grid.y_rho',1,2)); %diff(v,1,2)./repmat(diff(grid.y_v',1,2),[1 1 s(3) s(4)]);
    tz    = bsxfun(@rdivide,diff(lambda,1,3),permute(diff(grid.z_r,1,1),[3 2 1])); %diff(v,1,3)./repmat(permute(diff(grid.z_v,1,1),[3 2 1]),[1 1 1 s(4)]);
    
    % PV calculated at interior rho points
                                % f + vx - uy                      (rho)_z
    pv(:,:,:,tstart:tend) = double((avgx(avgz(bsxfun(@plus,avgy(vx - uy),f)))  .*  (tz(2:end-1,2:end-1,:,:)) ...
                   - avgy(vz(2:end-1,:,:,:).*avgz(ty(2:end-1,:,:,:))) ... % vz * (rho)_y
                   + avgx(uz(:,2:end-1,:,:).*avgz(tx(:,2:end-1,:,:))))./avgz(lambda(2:end-1,2:end-1,:,:))); % uz*(rho)_x
                               
end

% Write to netcdf

function [um] = avgy(um)
    um = (um(:,1:end-1,:,:)+um(:,2:end,:,:))/2;

function [um] = avgx(um)
    um = (um(1:end-1,:,:,:)+um(2:end,:,:,:))/2;

function [um] = avgz(um)
    um = (um(:,:,1:end-1,:)+um(:,:,2:end,:))/2;

    %% old code
    
%     pv1    = avgx(avgz(bsxfun(@plus,avgy(vx - uy),f)))  .*  (tz(2:end-1,2:end-1,:,:));
%     pv2    = (-1)*;
%     pv3    = uz.*avgz(tx);
    %pv = double((pv1 + avgy(pv2(2:end-1,:,:,:)) + avgx(pv3(:,2:end-1,:,:)))./avgz(lambda(2:end-1,2:end-1,:,:))); 