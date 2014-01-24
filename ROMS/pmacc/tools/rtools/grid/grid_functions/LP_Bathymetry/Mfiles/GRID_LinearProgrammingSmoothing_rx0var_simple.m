function NewBathy=GRID_LinearProgrammingSmoothing_rx0var_simple(...
    MSK, Hobs, rx0matrix)
% GRID_LinearProgrammingSmoothing_rx0_simple(MSK, Hobs, r, SignConst, AmpConst)
%
% ---MSK(eta_rho,xi_rho) is the mask of the grid
%      1 for sea
%      0 for land
% ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
% ---rx0matrix is the target roughness factor matrix
tol=0.0001;
K=find(MSK == 1);
if (min(Hobs(K)) < tol)
  disp('The bathymetry should always be positive');
  error('Please correct');
end;

[eta_rho, xi_rho]=size(Hobs);
SignConst=zeros(eta_rho, xi_rho);
AmpConst=10000*ones(eta_rho, xi_rho);

NewBathy=GRID_LinearProgrammingSmoothing_rx0(...
    MSK, Hobs, rx0matrix, SignConst, AmpConst);
    
