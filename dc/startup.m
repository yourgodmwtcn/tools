fprintf('\n startup.m entered \n');

% update tookboxcache
rehash toolboxcache


%% Make good figures

% figure properties
%set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor','w')
set(0,'DefaultFigureRenderer','painters')
set(0,'DefaultFigurePaperPositionMode', 'auto');
set(0,'DefaultTextFontName', 'AvantGarde');
set(0,'DefaultTextColor','k')
set(0,'DefaultTextFontSize',12)

% axes
set(0,'DefaultAxesFontName','AvantGarde')
set(0,'DefaultAxesTickLength'  , [.01 .01]);
set(0,'DefaultAxesLineWidth'  , 1);
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultAxesBox','on')
set(0,'DefaultAxesTickDir','in')
set(0,'DefaultAxesXMinorTick','on')
set(0,'DefaultAxesYMinorTick','on')
set(0,'DefaultAxesZMinorTick','on')
set(0,'DefaultAxesXColor',[.3 .3 .3])
set(0,'DefaultAxesYColor',[.3 .3 .3])
set(0,'DefaultAxesZColor',[.3 .3 .3])
set(0,'DefaultAxesLineWidth',1)

