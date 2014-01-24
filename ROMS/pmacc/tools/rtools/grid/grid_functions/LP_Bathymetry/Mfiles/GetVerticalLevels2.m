function [Z_r, Z_w]=GetVerticalLevels2(DEP_rho, MSK_rho, ARVD)
% [Z_r, Z_w]=GetVerticalLevels2(DEP_rho, MSK_rho, ARVD)
%
% Z_r, Z_w are vertical Z level (negative values)
% DEP_rho is depth at rho points
% MSK_rho is mask at rho points
% (it works also for u points, vpoints, psi-points)
%
% ARVD is a record containing the description of vertical 
% stratification. For example:
% ARVD.Vtransform=2;
% ARVD.Vstretching=1;
% ARVD.ThetaS=4;
% ARVD.ThetaB=0.35;
% ARVD.N=30;
% ARVD.hc=10;
%
% in output we should have 
% Z_r and Z_w negative
% Z_r interlaced inside Z_w
% Z_w(1,   :, :) is the bottom
% Z_w(N+1, :, :) is the surface
% if MSK_rho and DEP_rho are of size (eta_rho, xi_rho)
% then Z_r, Z_w are of size (N, eta_rho, xi_rho),
% (N+1, eta_rho, xi_rho)

[eta_rho, xi_rho]=size(DEP_rho);
if (isfield(ARVD, 'Zcoordinate') == 1)
  ListZ_r=ARVD.ListZ_r;
  ListZ_w=ARVD.ListZ_w;
  N=ARVD.N;
  Z_r=NaN*ones(N, eta_rho, xi_rho);
  Z_w=NaN*ones(N+1, eta_rho, xi_rho);
  for iEta=1:eta_rho
    for iXi=1:xi_rho
      if (MSK_rho(iEta, iXi) == 1)
	for i=1:N
	  if (-DEP_rho(iEta, iXi) < ListZ_r(i,1))
	    Z_r(i, iEta, iXi)=ListZ_r(i,1);
	    % if we have one wet r level, then we need wet w level
            % as well. It was a source of painful bugs.
	    Z_w(i, iEta, iXi)=ListZ_w(i,1);
	    Z_w(i+1, iEta, iXi)=ListZ_w(i+1,1);
	  end;
	end;
%	for i=1:N+1
%	  if (-DEP_rho(iEta, iXi) < ListZ_w(i,1))
%	    Z_w(i, iEta, iXi)=ListZ_w(i,1);
%	  end;
%	end;
      end;
    end;
  end;
elseif (isfield(ARVD, 'AREGcoordinate') == 1)
  ListZ_r=ARVD.ListZ_r;
  ListZ_w=ARVD.ListZ_w;
  N=ARVD.N;
  Z_r=zeros(N, eta_rho, xi_rho);
  Z_w=zeros(N+1, eta_rho, xi_rho);
  for iEta=1:eta_rho
    for iXi=1:xi_rho
      if (MSK_rho(iEta, iXi) == 1)
	for i=1:N
	  Z_r(i, iEta, iXi)=ListZ_r(i,1)*DEP_rho(iEta, iXi);
	end;
	for i=1:N+1
	  Z_w(i, iEta, iXi)=ListZ_w(i,1)*DEP_rho(iEta, iXi);
	end;
      end;
    end;
  end;
else
  [Sc_w, Cs_w, Sc_r, Cs_r]=GRID_GetSc_Cs_V2(ARVD);
  N=ARVD.N;
  DEPwork=DEP_rho;
  K=find(MSK_rho == 0);
  DEPwork(K)=3;
  Z_r=zeros(N, eta_rho, xi_rho);
  Z_w=zeros(N+1, eta_rho, xi_rho);
  if (ARVD.Vtransform == 1)
    for i=1:N
      Z_r(i,:,:)=ARVD.hc*Sc_r(i,1) + (DEPwork - ARVD.hc)*Cs_r(i,1);
    end;
    for i=1:N+1
      Z_w(i,:,:)=ARVD.hc*Sc_w(i,1) + (DEPwork - ARVD.hc)*Cs_w(i,1);
    end;
  elseif (ARVD.Vtransform == 2)
    for i=1:N
      Zo=(ARVD.hc*Sc_r(i,1) + DEPwork*Cs_r(i,1))./(ARVD.hc+DEPwork);
      Z_r(i,:,:)=Zo.*DEPwork;
    end;
    for i=1:N+1
      Zo=(ARVD.hc*Sc_w(i,1) + DEPwork*Cs_w(i,1))./(ARVD.hc+DEPwork);
      Z_w(i,:,:)=Zo.*DEPwork;
    end;
  else
    disp('Vtransform wrongly assigned');
    error('Please correct');
  end;
end;
