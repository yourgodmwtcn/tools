% calculates Ertel PV
%       [pv] = roms_pv(file,timesteps)

function [pv] = roms_pv(file,timesteps)

    file = find_file(file);

    %[vars,atts,dims] = ncdfread(file);
    grid = roms_get_grid(file,file,0,1);
    
    if ~exist('timesteps','var') || isempty(timesteps)
        timesteps = [1 Inf];
    end
    
    if numel(timesteps) == 1
        timesteps(2) = 1;
    end
    
    f      = ncread(file,'f',[1 1],[Inf Inf]);
    u      = ncread(file,'u',[1 1 1 timesteps(1)],[Inf Inf Inf (timesteps(2)-timesteps(1))]);
    v      = ncread(file,'v',[1 1 1 timesteps(1)],[Inf Inf Inf (timesteps(2)-timesteps(1))]);
    rho    = ncread(file,'rho',[1 1 1 timesteps(1)],[Inf Inf Inf (timesteps(2)-timesteps(1))]); % theta
        
    s = size(rho);
    if length(s) == 3
        s(4) = 1;
    end
    
    % get pot. temp.
    fmat   = repmat(f,[1 1 s(3) s(4)]);
    %pres   = repmat(reshape(grid.z_r,[s(1) s(2) s(3)]),[1 1 1 s(4)]);
    %theta  = sw_ptmp(32*ones(size(vars.temp)),vars.temp,pres,zeros(size(vars.temp)));
    lambda = rho; % choose either theta or rho
    
    % calculate gradients
    vx    = diff(v,1,1)./repmat(diff(grid.x_v',1,1),[1 1 s(3) s(4)]);
    vy    = diff(v,1,2)./repmat(diff(grid.y_v',1,2),[1 1 s(3) s(4)]);
    vz    = diff(v,1,3)./repmat(permute(diff(grid.z_v,1,1),[3 2 1]),[1 1 1 s(4)]);
    
    ux    = diff(u,1,1)./repmat(diff(grid.x_u',1,1),[1 1 s(3) s(4)]);
    uy    = diff(u,1,2)./repmat(diff(grid.y_u',1,2),[1 1 s(3) s(4)]);
    uz    = diff(u,1,3)./repmat(permute(diff(grid.z_u,1,1),[3 2 1]),[1 1 1 s(4)]);
    
    tx    = diff(lambda,1,1)./repmat(diff(grid.x_rho',1,1),[1 1 s(3) s(4)]);
    ty    = diff(lambda,1,2)./repmat(diff(grid.y_rho',1,2),[1 1 s(3) s(4)]);
    tz    = diff(lambda,1,3)./repmat(permute(diff(grid.z_r,1,1),[3 2 1]),[1 1 1 s(4)]);
    
    %% calculate pv
    coeff1 = (fmat(1:end-1,1:end-1,:,:) + vx - uy);
    pv1    = coeff1(:,:,1:end-1,:).*tz(1:end-1,1:end-1,:,:);
    pv2    = (-1)*vz.*ty(:,:,1:end-1,:);
    pv3    = uz.*tx(:,:,1:end-1,:);
    pv     = double((pv1 + pv2(1:end-1,:,:,:) + pv3(:,1:end-1,:,:))./rho(1:end-1,1:end-1,1:end-1,:));
    