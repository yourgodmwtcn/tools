% calculates Ertel PV
%       [pv] = roms_pv(file,timesteps)

function [pv] = roms_pv(file,timesteps)

% parameters
lam = 'rho';
vinfo = ncinfo(fname,'u');
dim   = length(vinfo.Size); 
slab  = 40;

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
f    = ncread(file,'f',[1 1],[Inf Inf]);
fmat   = repmat(f,[1 1 s(3) s(4)]);

for i=0:iend-1
    
    [read_start,read_count] = roms_ncread_params(dim,i,iend,slab,tindices,dt);

    u      = ncread(file,'u',read_start,read_count,stride);
    v      = ncread(file,'v',read_start,read_count,stride);
    lambda = ncread(file,lam,read_start,read_count,stride); % theta

    s = size(rho);
    if length(s) == 3
        s(4) = 1;
    end

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
    

    coeff1 = (fmat(1:end-1,1:end-1,:,:) + vx - uy);
    pv1    = coeff1(:,:,1:end-1,:).*tz(1:end-1,1:end-1,:,:);
    pv2    = (-1)*vz.*ty(:,:,1:end-1,:);
    pv3    = uz.*tx(:,:,1:end-1,:);
    pv     = double((pv1 + pv2(1:end-1,:,:,:) + pv3(:,1:end-1,:,:))./rho(1:end-1,1:end-1,1:end-1,:));
end


%     ux    = diff(u,1,1)./repmat(diff(grid.x_u',1,1),[1 1 s(3) s(4)]);
%     uy    = diff(u,1,2)./repmat(diff(grid.y_u',1,2),[1 1 s(3) s(4)]);
%     uz    = diff(u,1,3)./repmat(permute(diff(grid.z_u,1,1),[3 2 1]),[1 1 1 s(4)]);
%     
%     tx    = diff(lambda,1,1)./repmat(diff(grid.x_rho',1,1),[1 1 s(3) s(4)]);
%     ty    = diff(lambda,1,2)./repmat(diff(grid.y_rho',1,2),[1 1 s(3) s(4)]);
%     tz    = diff(lambda,1,3)./repmat(permute(diff(grid.z_r,1,1),[3 2 1]),[1 1 1 s(4)]);
% get pot. temp.

%pres   = repmat(reshape(grid.z_r,[s(1) s(2) s(3)]),[1 1 1 s(4)]);
%theta  = sw_ptmp(32*ones(size(vars.temp)),vars.temp,pres,zeros(size(vars.temp)));
%lambda = rho; % choose either theta or rho

%% calculate pv
