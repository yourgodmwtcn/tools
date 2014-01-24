function add_dyes_to_clm(clm_file,NESW,numdye)
% add_dyes_to_clm(clm_file,numdye)
% where numdye is the total number of dyes desired
% add dyes to the climatology files. Note that dye needs to be in these
% clm files but can be zero everywhere therefore the values never need to
% be filled in
%
% 12/20/2011 Sarah N. Giddings created adapted from Nick Lederer's code
% can be run from add_dyes
% e.g. clm_file = '/Users/sarahgid/Documents/Research/PNWTOX/runs/ptx_highT_2005_1/Ocn/ocean_clm_1.nc';

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
ncid = netcdf.open(clm_file,'NC_WRITE');
%set the fill value to nothing so it will be zeros by default
netcdf.setFill(ncid,'NC_NOFILL');

%get the salt time id
salt_time_dim_id = netcdf.inqDimID(ncid,'salt_time');
%get the info about the salt time including dimensions to copy to dye time
[~, salt_time_len] = netcdf.inqDim(ncid, salt_time_dim_id);

%get the salt id
salt_id = netcdf.inqVarID(ncid,'salt');
%get the info about the salt including dimensions to copy to dye
[~, xtype, dimids, ~] = netcdf.inqVar(ncid,salt_id);

%put the netcdf file in define mode
netcdf.reDef(ncid);

%define the dye time id
dye_time_dim_id = netcdf.defDim(ncid,'dye_time',salt_time_len);
%define the dye time based on the salt time dimensions and type
dimids(dimids==salt_time_dim_id) = dye_time_dim_id;

%loop through each dye
for i=1:num_dye
    varname = sprintf('dye_%02d',i);    
    %define the variable based on the salt dimensions and type
    varid = netcdf.defVar(ncid,varname,xtype,dimids);
    %add the desired attributes
    netcdf.putAtt(ncid,varid,'long_name','dye concentration');
    netcdf.putAtt(ncid,varid,'units','kg/m3');
end

% dye_time variable definition
salt_time_id = netcdf.inqVarID(ncid,'salt_time');
[~, xtype, dimids, ~] = netcdf.inqVar(ncid,salt_time_id);
dimids(dimids==salt_time_dim_id) = dye_time_dim_id;
dye_time_id =  netcdf.defVar(ncid,'dye_time',xtype,dimids);
%add the desired attributes
netcdf.putAtt(ncid,dye_time_id,'long_name','dye concentration climatology time');
netcdf.putAtt(ncid,dye_time_id,'units','seconds');

%finish defining
netcdf.endDef(ncid);

% copy salt_time into dye_time
salt_time = netcdf.getVar(ncid,salt_time_id);
netcdf.putVar(ncid,dye_time_id,salt_time);

%loop through each boundary dye adding ones to the climatology file such
%that nudging to climatology will be sort of like nudging to tracers? Note
%that this is not a great solution, better would be to change
%ana_nudgcoef.h significantly OR upgrade to the new ROMS where boundary
%conditions can be set differently for each tracer in the .in file
for i = num_dye_rivers + 1:num_dye
    varname = sprintf('dye_%02d',i);
    dye_id = netcdf.inqVarID(ncid,varname);
    dye = netcdf.getVar(ncid,dye_id);
    dye(:) = 1;
    netcdf.putVar(ncid,dye_id,dye);
end
    
%close the netcdf file
netcdf.close(ncid);
tclm = toc;
disp(['Completed adding ' num2str(num_dye) ' dye variables to clm file '...
    clm_file ' in ' num2str(tclm) ' seconds'])


