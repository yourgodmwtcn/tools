function add_dyes_to_rivers(river_file,AllRiverNames,DyeRiverNames,numdye)
% add_dyes_to_rivers(river_file,riverRefs,numdye)
% where river_file is the rivers.nc file with rivers in the order listed
% AllRiverNames contains a string list of ALL river names in river_file
% DyeRiverNames contains a string list of river names which you want to add dye
% e.g. AllRiverNames =
% {'skagit','snohomish','stillaguamish','columbia','puyallup','duwamish','nisqually','deschutes','skagit_south','skokomish','duckabush','fraser','dosewallips','hammahamma','cedar','nooksack','samish'};
% e.g. DyeRiverNames = {'columbia','fraser','others'};
% numdye contains the total number of dyes to input
%
% 12/20/2011 Sarah N. Giddings created adapted from Nick Lederer's code
% and updated to make more general, again more general 1/3/2012

tic
%number of river dyes
num_dye_rivers = length(DyeRiverNames);
%number of boundary dyes
num_dye_brys = numdye - num_dye_rivers;

%open the netcdf file
ncid = netcdf.open(river_file,'NC_WRITE');
%set the fill value to nothing so it will be zeros by default
netcdf.setFill(ncid,'NC_NOFILL');

%get the salt id
salt_id = netcdf.inqVarID(ncid,'river_salt');
%get the info about the salt including dimensions to copy to dye
[~, xtype, dimids, ~] = netcdf.inqVar(ncid,salt_id);

%list of all rivers which gets reduced as individual river dyes are added
%thus leaving a list of all other rivers without dye
otherRivers = 1:length(AllRiverNames);
%loop through each river dye
for i=1:num_dye_rivers
    %turn on defining
    netcdf.reDef(ncid);
    varname = sprintf('river_dye_%02d',i);
    %define the variable based on the salt dimensions and type
    varid = netcdf.defVar(ncid,varname,xtype,dimids);
    %add the desired attributes
    netcdf.putAtt(ncid,varid,'long_name','river dye concentration');
    netcdf.putAtt(ncid,varid,'units','kg/m3');
    %finish defining
    netcdf.endDef(ncid);
    
    %get the dye (all zeros) and fill it in with 1 where appropriate
    dye = netcdf.getVar(ncid,varid);
    %loop through each river in AllRiverNames to ensure you are adding dye
    %to the appropriate river number
    for j = 1:length(AllRiverNames)
        if strcmp(DyeRiverNames(i),AllRiverNames(j))
            dye(j,:,:) = 1;
            %remove river indices that already have dye from otherRivers
            otherRivers(otherRivers==j) = [];
            disp(['adding dye ' num2str(j) ' to ' AllRiverNames(j)])
        end
    end
    if strcmp(DyeRiverNames(i),'others')
        % All other rivers
        dye(otherRivers,:,:) = 1;
        disp(['adding dye ' num2str(otherRivers) ' to remaining rivers'])
    end
    %put the variable in the netcdf file
    netcdf.putVar(ncid,varid,dye);
end

%loop through each boundary dye
for i=num_dye_rivers+1:numdye
    %turn on defining
    netcdf.reDef(ncid);
    varname = sprintf('river_dye_%02d',i);
    %define the variable based on the salt dimensions and type
    varid = netcdf.defVar(ncid,varname,xtype,dimids);
    %add the desired attributes
    netcdf.putAtt(ncid,varid,'long_name','boundary dye concentration');
    netcdf.putAtt(ncid,varid,'units','kg/m3');
    %finish defining
    netcdf.endDef(ncid);
end

%close the netcdf file
netcdf.close(ncid);    
tclm = toc;
disp(['Completed adding ' num2str(num_dye_rivers) ' river dye variables and ' ...
    num2str(num_dye_brys) ' boundary dye variables to river file '...
    river_file ' in ' num2str(tclm) ' seconds'])


    