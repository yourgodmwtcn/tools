function NewBathy=GRID_LinearProgrammingSmoothing_rx0_similarity(...
    MSK, Hobs, rx0max, SignConst, AmpConst, MaxSimilarity)
% GRID_LinearProgrammingSmoothing_rx0_similarity(...
%    MSK, Hobs, rx0max, SignConst, AmpConst, MaxSimilarity)
%
% ---MSK(eta_rho,xi_rho) is the mask of the grd
%      1 for sea
%      0 for land
% ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
% ---rx0max is the target rx0 roughness factor
% ---SignConst(eta_rho,xi_rho) matrix of 0, +1, -1
%      +1  only bathymetry increase are allowed.
%      -1  only bathymetry decrease are allowed.
%      0   increase and decrease are allowed.
%      (put 0 if you are indifferent)
% ---AmpConst(eta_rho,xi_rho) matrix of reals.
%      coefficient alpha such that the new bathymetry should
%      satisfy to  |h^{new} - h^{raw}| <= alpha h^{old}
%      (put 10000 if you are indifferent)
% ---MaxSimilarity is the maximal expansion factor on every
%    cell
[eta_rho, xi_rho]=size(MSK);
rx0matrix=rx0max*ones(eta_rho, xi_rho);
NewBathy=GRID_LinearProgrammingSmoothing_rx0var_similarity(...
    MSK, Hobs, rx0matrix, SignConst, AmpConst, MaxSimilarity);
