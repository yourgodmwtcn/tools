function varlist = roms_varlist(category,model)
% get cell array of subsets of roms variable names
%
% varlist = roms_varlist(category,model)
%
% categories are: physics, physics2d, physics3d, mixing3d, s-param, 
% s-coord, grid, fennel, dom, biodiags, bulkflux

% $Id: roms_varlist.m 407 2012-02-16 15:44:10Z wilkin $

if nargin == 0
  type(mfilename)
end

if nargin < 2
  model = 'roms';
end

switch model
  
  case 'roms'
    switch category
      case 'time'
        varlist = 'ocean_time';
      case 'physics'
        varlist = {'temp','salt','u','v','zeta','ubar','vbar'};
      case 'physics2d'
        varlist = {'zeta','ubar','vbar'};
      case 'physics3d'
        varlist = {'temp','salt','u','v'};
      case 'mixing3d'
        varlist = {'AKv','AKt','AKs'};
      case 's-param'
        varlist = {'theta_s','theta_b','Tcline','hc','Vtransform','Vstretching'};
      case 's-coord'
        varlist = {'s_rho','s_w','Cs_r','Cs_w'};
      case 'grid'
        varlist = {'h','f','pm','pn','angle','lon_rho','lat_rho',...
          'lon_u','lat_u','lon_v','lat_v','lon_psi','lat_psi',...
          'mask_rho','mask_u','mask_v','mask_psi'};
      case 'fennel'
        varlist = {'NO3','NH4','chlorophyll','phytoplankton','zooplankton',...
          'LdetritusN','SdetritusN','TIC','alkalinity','LdetritusC',...
          'SdetritusC','oxygen'};
      case {'dom','usecos'}
        varlist = {'semilabileDON','refractoryDON',...
          'semilabileDOC','refractoryDOC'};
      case 'biodiags'
        varlist = {'denitrification','CO2_airsea','pCO2',...
          'P_Production','NO3_uptake',...
          'nitrogen_buried','carbon_bottom','carbon_buried','C_excess_uptake'};
      case 'bulkflux'
        varlist = {'Uwind','Vwind','Pair','Tair','Qair','swrad','lwrad','lwrad_down','rain'};
    end
    
  case 'ncom'
    switch category
      case 'physics'
        varlist = {'water_temp','salinity','water_u','water_v','surf_el',[],[]};
      case 'physics2d'
        varlist = {'surf_el',[],[]};
      case 'physics3d'
        varlist = {'water_temp','salinity','water_u','water_v'};
    end
    
end
