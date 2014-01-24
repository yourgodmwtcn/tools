function add_dyes_to_bry(bry_file,NESW,numdye)
% add_dyes_to_bry(bry_file,NESW,numdye)
% where NESW is a vector of 4 values indicating whether to add dye to the
% northern, eastern, southern, and western boundaries respectively
% and where numdye is the total number of dyes desired
%
% 12/20/2011 Sarah N. Giddings created adapted from Nick Lederer's code
% can be run from add_dyes
% e.g. bry_file = '/Users/sarahgid/Documents/Research/PNWTOX/runs/ptx_highT_2005_1/Ocn/ocean_bry_1.nc';
% NESW = [0 0 1 1]; will put dye in the outhern and western boundaries
% only

if nargin < 2
    NESW = [0 0 0 0]; %default is no boundary dye
elseif nargin == 3
    num_dye = numdye;
else
    num_dye = 5;
end
%number of boundary dyes
num_dye_brys = nansum(NESW);
%number of river dyes
num_dye_rivers = num_dye - num_dye_brys;

tic
%open the netcdf file
ncid = netcdf.open(bry_file,'NC_WRITE');
%set the fill value to nothing so it will be zeros by default
netcdf.setFill(ncid,'NC_NOFILL');

%get the salt time dimension information
salt_time_dim_id = netcdf.inqDimID(ncid,'salt_time');
[~, salt_time_len] = netcdf.inqDim(ncid, salt_time_dim_id);
netcdf.reDef(ncid);

%define the dye time dimension
dye_time_dim_id = netcdf.defDim(ncid,'dye_time',salt_time_len);

%set the bry values at each boundary to 0 throughout depth, time, and along
%the full length of the boundary
% NORTH
salt_north_id = netcdf.inqVarID(ncid,'salt_north');
[~, xtype, dimids, ~] = netcdf.inqVar(ncid,salt_north_id);
dimids(dimids==salt_time_dim_id) = dye_time_dim_id;

for i=1:num_dye
    varname = sprintf('dye_north_%02d',i);
    netcdf.defVar(ncid,varname,xtype,dimids);
end

% SOUTH
salt_south_id = netcdf.inqVarID(ncid,'salt_south');
[~, xtype, dimids, ~] = netcdf.inqVar(ncid,salt_south_id);
dimids(dimids==salt_time_dim_id) = dye_time_dim_id;

for i=1:num_dye
    varname = sprintf('dye_south_%02d',i);
    netcdf.defVar(ncid,varname,xtype,dimids);
end

% EAST
salt_east_id = netcdf.inqVarID(ncid,'salt_east');
[~, xtype, dimids, ~] = netcdf.inqVar(ncid,salt_east_id);
dimids(dimids==salt_time_dim_id) = dye_time_dim_id;

for i=1:num_dye
    varname = sprintf('dye_east_%02d',i);
    netcdf.defVar(ncid,varname,xtype,dimids);
end

% WEST
salt_west_id = netcdf.inqVarID(ncid,'salt_west');
[~, xtype, dimids, ~] = netcdf.inqVar(ncid,salt_west_id);
dimids(dimids==salt_time_dim_id) = dye_time_dim_id;

for i=1:num_dye
    varname = sprintf('dye_west_%02d',i);
    netcdf.defVar(ncid,varname,xtype,dimids);
end

% dye_time variable definition
salt_time_id = netcdf.inqVarID(ncid,'salt_time');
[~, xtype, dimids, ~] = netcdf.inqVar(ncid,salt_time_id);
dimids(dimids==salt_time_dim_id) = dye_time_dim_id;
dye_time_id =  netcdf.defVar(ncid,'dye_time',xtype,dimids);

netcdf.endDef(ncid);

% copy salt_time into dye_time
salt_time = netcdf.getVar(ncid,salt_time_id);
netcdf.putVar(ncid,dye_time_id,salt_time);

%get the dye (all zeros) and fill it in with 1 where appropriate
%start with the first dye number after the river dyes and then increment up
%from there
i = num_dye_rivers + 1;
if NESW(1) == 1
    % NORTH
    varname = sprintf('dye_north_%02d',i);
    dye_id = netcdf.inqVarID(ncid,varname);
    dye = netcdf.getVar(ncid,dye_id);
    dye(:,:,:) = 1;
    netcdf.putVar(ncid,dye_id,dye);
    disp(['added dye ' varname ' to northern boundary'])
    i = i+1;
end
if NESW(2) == 1
    % EAST
    varname = sprintf('dye_east_%02d',i);
    dye_id = netcdf.inqVarID(ncid,varname);
    dye = netcdf.getVar(ncid,dye_id);
    dye(:,:,:) = 1;
    netcdf.putVar(ncid,dye_id,dye);
    disp(['added dye ' varname ' to eastern boundary'])
    i = i+1;
end
if NESW(3) == 1
    % SOUTH
    varname = sprintf('dye_south_%02d',i);
    dye_id = netcdf.inqVarID(ncid,varname);
    dye = netcdf.getVar(ncid,dye_id);
    dye(:,:,:) = 1;
    netcdf.putVar(ncid,dye_id,dye);
    disp(['added dye ' varname ' to southern boundary'])
    i = i+1;
end
if NESW(4) == 1
    % WEST
    varname = sprintf('dye_west_%02d',i);
    dye_id = netcdf.inqVarID(ncid,varname);
    dye = netcdf.getVar(ncid,dye_id);
    dye(:,:,:) = 1;
    netcdf.putVar(ncid,dye_id,dye);
    disp(['added dye ' varname ' to western boundary'])
    i = i+1;
end

%close the netcdf file
netcdf.close(ncid);
tclm = toc;
disp(['Completed adding ' num2str(num_dye_rivers) ' river dye variables and ' ...
    num2str(num_dye_brys) ' boundary dye variables to bry file '...
    bry_file ' in ' num2str(tclm) ' seconds'])

