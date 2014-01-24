FileSave='GRID_ADRIA02';
load(FileSave, 'mreza');
MSK_rho=mreza.MSK_rho;
DEP_rho=mreza.SampledBathy;
[eta_rho, xi_rho]=size(DEP_rho);


rx0max=0.1;
disp('Using LP with varying rx0');
KseaRho=find(MSK_rho == 1);
nbWet=size(KseaRho,1);

avgDep=sum(mreza.SampledBathy(KseaRho))/nbWet;
TheMult=(avgDep + DEP_rho)/avgDep;
K=find(TheMult > 2);
TheMult(K)=2;
rx0matrix=rx0max*TheMult;



NewBathy=GRID_LinProgHeuristic_rx0var(...
    MSK_rho, DEP_rho, rx0matrix);
