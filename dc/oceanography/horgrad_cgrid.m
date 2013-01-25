% calculates horizontal gradient of variable 'var' along axis 'ax'
% does central difference
%   grid - has grid.xmat,grid.ymat,grid.zmat
%   var  - variable
%   ax - 1,2 for x,y

function [horgrad] = horgrad_cgrid(grid,var,ax)

% #1 - loop through every point, linearly interpolate neighbouring
% profiles to same depth and calculate derivative
s = size(var);

if ax == 1
    axmat = grid.xmat;
else
    axmat = grid.ymat;
end

horgrad = nan(s);

for ii = 2:s(1)-1
    for jj = 2:s(2)-2
        % appropriate indices
        if ax == 1
            im1 = ii-1;
            ip1 = ii+1;
            jm1 = jj;
            jp1 = jj;
        else
            im1 = ii;
            ip1 = ii;
            jm1 = jj-1;
            jp1 = jj+1;
        end
        
        % figure out appropriate z grids
        zz   = grid.zmat(ii,jj,:);
        zzm1 = grid.zmat(im1,jm1,:);
        zzp1 = grid.zmat(ip1,jp1,:);
        
        for kk = 1:s(3)
            % var @ minus 1 point interpolated
            vm1i = interp1(zzm1,squeeze(var(im1,jm1,:)),zz);
            vp1i = interp1(zzp1,squeeze(var(ip1,jp1,:)),zz);
            
            dx   = (axmat(ip1,jp1,kk) - axmat(im1,jm1,kk))/2;
           
            horgrad(ii,jj,kk) = (vp1i(kk)-vm1i(kk))./dx;
        end
    end
end

% #2 - interpolate to uniform grid, differentiate and then interpolate back
% wont work because i'll need a lot of points to preserve resolution in
% shallow water
% zi = linspace(max(grid.zmat(:)),min(grid.zmat(:)),100);
% [ximat,yimat,zimat] = ndgrid(grid.xmat(:,1,1),grid.ymat(1,:,1),zi);
% 
% vari = interpn(grid.xmat,grid.ymat,grid.zmat,var,ximat,yimat,zimat);