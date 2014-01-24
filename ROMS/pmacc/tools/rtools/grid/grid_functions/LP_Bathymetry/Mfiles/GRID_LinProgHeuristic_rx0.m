function TheNewBathy=GRID_LinProgHeuristic_rx0(...
    MSK_rho, DEP_rho, rx0max)


% We are doing smoothing with respect to rx0
% The point of this version is that we are using an heuristic
% that speed up things considerably for the run.

[eta_rho, xi_rho]=size(DEP_rho);

rx0matrix=rx0max*ones(eta_rho, xi_rho);
TheNewBathy=GRID_LinProgHeuristic_rx0var(...
    MSK_rho, DEP_rho, rx0matrix);

