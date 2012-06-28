function [pm,pn] = Z_pmpn(x_u,y_v);

% 10/17/2008  Parker MacCready
% This returns the arrays of metrics, pm and pn.

dx = diff(x_u,1,2);
dx_rho = [dx(:,1), dx, dx(:,end)];
dy = diff(y_v,1,1);
dy_rho = [dy(1,:); dy; dy(end,:)];
pm = 1./dx_rho;
pn = 1./dy_rho;
