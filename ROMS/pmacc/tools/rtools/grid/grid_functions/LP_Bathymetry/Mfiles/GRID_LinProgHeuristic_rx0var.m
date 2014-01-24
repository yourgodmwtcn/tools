function TheNewBathy=GRID_LinProgHeuristic_rx0var(...
    MSK_rho, DEP_rho, rx0matrix)

% the points that need to be modified in the 
% version by 



% paperwork
KseaRho=find(MSK_rho == 1);
nbWet=size(KseaRho,1);
[eta_rho, xi_rho]=size(DEP_rho);
disp(['eta_rho=' num2str(eta_rho) ' xi_rho=' num2str(xi_rho) ...
      ' nbWet=' num2str(nbWet)]);
ETAmat=zeros(eta_rho, xi_rho);
XImat=zeros(eta_rho, xi_rho);
for iEta=1:eta_rho
  for iXi=1:xi_rho
    ETAmat(iEta, iXi)=iEta;
    XImat(iEta, iXi)=iXi;
  end;
end;

% definition of problem
RMat=GRID_RoughnessMatrix(DEP_rho, MSK_rho);
DiffRX0=rx0matrix - RMat;
TheMin=min(DiffRX0(KseaRho));
disp([' input: min(rx0matrix - RMat)=' num2str(TheMin)]);

Kbad=find(RMat > rx0matrix);
MSKbad=zeros(eta_rho, xi_rho);
MSKbad(Kbad)=1;
nbKbad=size(Kbad,1);
disp(['found nbKbad=' num2str(nbKbad) ' bad points']);




Kdist=5;
disp(['Determining neighborhood graph for Kdist=' num2str(Kdist)]);
ListIdx=zeros(eta_rho, xi_rho);
ListIdx(Kbad)=1:nbKbad;

ETA_K=ETAmat(Kbad);
XI_K=XImat(Kbad);
ListEdges=zeros(0,2);
nbEdge=0;
for iK=1:nbKbad
  iVert=Kbad(iK,1);
  iEta=ETAmat(iVert);
  iXi=XImat(iVert);
  ListNeigh=SUB_StepNeighborhood(MSK_rho, iEta, iXi, 2*Kdist+1);
  nbNeigh=size(ListNeigh,1);
  for iNeigh=1:nbNeigh
    iEtaN=ListNeigh(iNeigh,1);
    iXiN=ListNeigh(iNeigh,2);
    if (MSKbad(iEtaN, iXiN) == 1)
      idx=ListIdx(iEtaN, iXiN);
      if (idx > iK)
	nbEdge=nbEdge+1;
	ListEdges(nbEdge,1)=iK;
	ListEdges(nbEdge,2)=idx;
      end;
    end;
  end;
end;

disp(['Determining connected components of the graph']);
ListVertexStatus=GRAPH_ConnectedComponent(...
    ListEdges, nbKbad);
nbColor=max(ListVertexStatus);
for iColor=1:nbColor
  K=find(ListVertexStatus == iColor);
  nbK=size(K,1);
  disp(['iColor=' num2str(iColor) '  nbK=' num2str(nbK)]);
end;
TheNewBathy=DEP_rho;
for iColor=1:nbColor
  disp('---------------------------------------------------------------');
  MSKcolor=zeros(eta_rho, xi_rho);
  K=find(ListVertexStatus == iColor);
  nbK=size(K,1);
  idx=0;
  for iVertex=1:nbKbad
    if (ListVertexStatus(iVertex,1) == iColor)
      idx=idx+1;
      iVert=Kbad(iVertex,1);
      iEta=ETAmat(iVert);
      iXi=XImat(iVert);
      MSKcolor(iEta, iXi)=1;
      ListNeigh=SUB_StepNeighborhood(MSK_rho, iEta, iXi, Kdist);
      nbNeigh=size(ListNeigh,1);
      for iNeigh=1:nbNeigh
	iEtaN=ListNeigh(iNeigh,1);
	iXiN=ListNeigh(iNeigh,2);
	MSKcolor(iEtaN, iXiN)=1;
      end;
    end;
  end;
  Kcolor=find(MSKcolor == 1);
  nbKcolor=size(Kcolor,1);
  disp(['iColor=' num2str(iColor) '  nbK=' num2str(nbK) ...
	' nbKcolor=' num2str(nbKcolor)]);
  Hobs=zeros(eta_rho, xi_rho);
  Hobs(Kcolor)=DEP_rho(Kcolor);
  NewBathy=GRID_LinearProgrammingSmoothing_rx0var_simple(...
      MSKcolor, Hobs, rx0matrix);
  TheNewBathy(Kcolor)=NewBathy(Kcolor);
end;
disp('Final obtained bathymetry');
RMat=GRID_RoughnessMatrix(TheNewBathy, MSK_rho);
DiffRX0=rx0matrix - RMat;
TheMin=min(DiffRX0(KseaRho));
disp(['output: min(rx0matrix - RMat)=' num2str(TheMin)]);
