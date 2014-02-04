% (R)oms (N)umerical (T)oolbox
%
% FUNCTION grdinfo = rnt_gridinfo(gridid)
%
% Loads the grid configuration for gridid
% To add new grid please edit this file.
% just copy an existing one and modify for
% your needs. It is simple.
%
% If you editing this file after using
% the Grid-pak scripts use the content
% of variable "nameit" for gridid.
%
% Example: CalCOFI application
%
%    grdinfo = rnt_gridinfo('calc')
%
% RNT - E. Di Lorenzo (edl@gatech.edu)

function gridindo=rnt_gridinfo(gridid)
  
% initialize to defaults
  gridindo.id      = gridid;
  gridindo.name    = '';
  gridindo.grdfile = '';
  gridindo.N       = 20;
  gridindo.thetas  = 5;
  gridindo.thetab  = 0.4;
  gridindo.tcline  = 200;
  gridindo.cstfile = which('rgrd_WorldCstLinePacific.mat');
  
  if exist(gridid)== 2
     file = textread(gridid,'%s','delimiter','\n','whitespace','');
     for i=1:length(file)
       eval(file{i});
     end
%     feval(gridid);
%     load(gridid);
     return
  end
  
  
  switch gridid


  case 'wrd'
    gridindo.id      = gridid;
    gridindo.name    = 'World';
    gridindo.grdfile = which('wrd-grid.nc');
    gridindo.N       = 10;
    gridindo.thetas  = 3;
    gridindo.thetab  = 0.4;
    gridindo.tcline  = 200;
    gridindo.tcline  = 50;

  case 'ias20'
    gridindo.id      = gridid;
    gridindo.name    = 'Intra-america-sea 20km';
    gridindo.grdfile = which('ias20_grid.nc');
    gridindo.N       = 30;
    gridindo.thetas  = 5;
    gridindo.thetab  = 0.4;
    gridindo.hc  = 20.;
    gridindo.cstfile ='ias_coast.mat';
    
  case 'northsea'
    gridindo.id      = gridid;
    gridindo.name    = 'North Sea';
    gridindo.grdfile = which('northsea-grid.nc');
    gridindo.N       = 25;
    gridindo.thetas  = 3;
    gridindo.thetab  = 0.4;
    gridindo.cstfile = which('rgrd_CoastlineWorld.mat');
            
  case 'china25'
    gridindo.id      = gridid;
    gridindo.name    = 'East China Sea 25 km';
    gridindo.grdfile = which('china-grid.nc');
    gridindo.N       = 20;
    gridindo.thetas  = 5;
    gridindo.thetab  = 0.4;

  case 'pac25'
    gridindo.id      = gridid;
    gridindo.name    = 'Pacific';
    gridindo.grdfile = which('pac25-grid.nc');
    gridindo.N       = 40;
    gridindo.thetas  = 6;
    gridindo.hc  = 50;
    gridindo.thetab  = 0.0;

  case 'indian'
    gridindo.id      = gridid;
    gridindo.name    = 'Indian Ocean';
    gridindo.grdfile = which('indian-grid.nc');
    gridindo.N       = 20;
    gridindo.thetas  = 6;
    gridindo.thetab  = 0.0;

case 'indo-pac'
    gridindo.id      = gridid;
    gridindo.name    = 'Indo-pacific';
    gridindo.grdfile = which('indo-pac-grid.nc');
    gridindo.N       = 20;
    gridindo.thetas  = 6;
    gridindo.thetab  = 0.0;

  case 'nepd'
    gridindo.id      = gridid;
    gridindo.name    = 'NEPD-GOA-CCS grid';
    gridindo.grdfile = which('nepd-grid.nc');
    gridindo.N       = 30;
    gridindo.thetas  = 5;
    gridindo.thetab  = 0.4;

  case 'scb'
    gridindo.id      = gridid;
    gridindo.name    = 'SCB 7 km';
    gridindo.grdfile = which('scb-grid.nc');
    gridindo.N       = 20;
    gridindo.thetas  = 6.5;
    gridindo.thetab  = 0;
    gridindo.hc  = 155.6225;
    gridindo.cstfile = which('rgrd_CCS_CstLine.mat');

  case 'nep10'
    gridindo.id      = gridid;
    gridindo.name    = 'NEP 10 km new';
    gridindo.grdfile = which('nep10-grid.nc');
    gridindo.N       = 42;
    gridindo.thetas  = 5;
    gridindo.thetab  = 0.4;
    gridinfo.Tcline  = 50;
    gridinfo.hc      = 30;

  case 'ccs10'
    gridindo.id      = gridid;
    gridindo.name    = 'CCS extract of NEP 10 km new';
    gridindo.grdfile = which('ccs10-grid.nc');
    gridindo.N       = 42;
    gridindo.thetas  = 5;
    gridindo.thetab  = 0.4;
    gridinfo.Tcline  = 50;
    gridinfo.hc      = 30;


  case 'nepd-ccs'
    gridindo.id      = gridid;
    gridindo.name    = 'NEPD-CCS grid';
    gridindo.grdfile = which('nepd-ccs-grid.nc');
    gridindo.N       = 30;
    gridindo.thetas  = 5;
    gridindo.thetab  = 0.4;

  case 'nepd2x'
    gridindo.id      = gridid;
    gridindo.name    = 'CCS grid';
    gridindo.grdfile = which('nepd2x-grid.nc');
    gridindo.N       = 30;
    gridindo.thetas  = 5;
    gridindo.thetab  = 0.4;    
    
  otherwise
    gridindo.id      = gridid;
    gridindo.name    = 'null';
    gridindo.grdfile = '/dev/null';
    gridindo.N       = 0;
    gridindo.thetas  = 0;
    gridindo.thetab  = 0;
    gridindo.tcline  = 0;
    disp([' RNT_GRIDINFO - ',gridid,' not configured']);
  end

% /neo/ROMS_Tutorial/matlib/rnt/rnt_listgrids.m
