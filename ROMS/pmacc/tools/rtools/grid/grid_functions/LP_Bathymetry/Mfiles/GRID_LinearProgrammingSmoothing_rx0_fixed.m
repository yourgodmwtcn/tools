function NewBathy=GRID_LinearProgrammingSmoothing_rx0_fixed(...
    MSK, Hobs, PRS, rx0max)
% this is a program for doing linear programming optimization
% of the bathymetry with some point being fixed.
% GRID_LinearProgrammingSmoothing_fixed(MSK, Hobs, PRS, r)
%
% ---MSK(eta_rho,xi_rho) is the mask of the grd
%      1 for sea
%      0 for land
% ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
% ---PRS(eta_rho,xi_rho) is the array for preservation status.
%      1 for preserving depth
%      0 for not preserving depth
% ---rx0max is the target rx0 roughness factor

[eta_rho,xi_rho]=size(MSK);
rx0matrix=rx0max*ones(eta_rho,xi_rho);
NewBathy=GRID_LinearProgrammingSmoothing_rx0var_fixed(...
    MSK, Hobs, PRS, rx0matrix);
