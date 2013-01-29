% calculates horizontal gradient of variable 'var' along axis 'ax' using
% co-ordinate transformation as in WikiROMS
% does _forward_ difference 
%   rgrid - from roms_get_grid.m - IMPORTANT - expects [z,y,x]
%         - needs rgrid.z_w, rgrid.s_w
%   vgrid - has vgrid.xmat,vgrid.ymat,vgrid.zmat for var - expects [x,y,z]
%           and also appropriate vgrid.s vector (S.s_w or S.s_rho)
%   var  - variable
%   ax1 - 1,2 for x,y

function [horgrad] = horgrad_cgrid(rgrid,vgrid,var,ax1)

%[s,C] = stretching(rgrid.Vstretching,rgrid.theta_s,rgrid.theta_b,rgrid.hc,rgrid.N,0,0);

switch ax1
    case 1
        axmat = vgrid.xmat;
    case 2 
        axmat = vgrid.ymat;
end

% Hz = dz/d?
% We compute Hz discretely as ?z/?? since this leads to the vertical sum of Hz
% being exactly the total water depth D. (from WikiROMS / Manual)
Hz = permute(bsxfun(@rdivide,diff(rgrid.z_w,1,1),diff(rgrid.s_w')),[3 2 1]);

% (dz/dx)_?
dzdx_s = diff(vgrid.zmat,1,ax1)./diff(axmat,1,ax1);

% (df/dx)_?; x = ax1
dfdx_s = diff(var,1,ax1)./diff(axmat,1,ax1);

% df/d?
dfds = nan(size(var));
dfds(:,:,2:end-1) = avg1(bsxfun(@rdivide,diff(var,1,3),permute(diff(vgrid.s'),[3 2 1])),3);
% pad on bottom and surface values.
dfds(:,:,1) = dfds(:,:,2);
dfds(:,:,end) = dfds(:,:,end-1); 

% chain rule power!
horgrad = dfdx_s - avg1(1./Hz,1) .* dzdx_s .* avg1(dfds,1);
