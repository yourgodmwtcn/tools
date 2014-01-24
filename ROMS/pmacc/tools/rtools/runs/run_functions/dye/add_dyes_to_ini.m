function add_dyes_to_ini(ini_file,numdye)
% add_dyes_to_ini(ocean_dir,numdye)
% where numdye is the total number of dyes desired
% add dyes to the initial conditions file. Note that dye needs to be in
% this file but can be zero everywhere therefore the values never need to
% be filled in
%
% 12/20/2011 Sarah N. Giddings created adapted from Nick Lederer's code
% can be run from add_dyes
% e.g. ini_file = '/Users/sarahgid/Documents/Research/PNWTOX/runs/ptx_highT_2005_1/Ocn/ocean_ini_1.nc';

if nargin == 2
    num_dye = numdye;
else
    num_dye = 5;
end

tic
%open the netcdf file
ncid = netcdf.open(ini_file,'NC_WRITE');
%set the fill value to nothing so it will be zeros by default
netcdf.setFill(ncid,'NC_NOFILL');
%get the salt id
salt_id = netcdf.inqVarID(ncid,'salt');
%get the info about the salt including dimensions to copy to dye
[~, xtype, dimids, ~] = netcdf.inqVar(ncid,salt_id);

%put the netcdf file in define mode
netcdf.reDef(ncid);
%loop through each dye
for i=1:num_dye
    varname = sprintf('dye_%02d',i);
    %define the variable based on the salt dimensions and type
    varid = netcdf.defVar(ncid,varname,xtype,dimids);
    %add the desired attributes
    netcdf.putAtt(ncid,varid,'long_name','dye concentration');
    netcdf.putAtt(ncid,varid,'units','kg/m3');
end
%finish defining
netcdf.endDef(ncid);
%close the netcdf file
netcdf.close(ncid);
tclm = toc;
disp(['Completed adding ' num2str(num_dye) ' dye variables to ini file '...
    ini_file ' in ' num2str(tclm) ' seconds'])



