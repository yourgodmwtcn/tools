function[vars varatts dims] = ncdfread(fname)

% read a netcdf file in the new MATLAB 7.11 PITA way
%       [vars varatts dims] = ncdfread(fname)
% output:
%       vars: the variables and their values
%       varatts: the variable attributes and their values
%       dims: dimensions of all variables

% currently set up to read, not write a netcdf file

% CHANGELOG Deepak Cherian
% Added wildcard matching of filename                      16th Jan 2012
% Modified to fill in NaN's, scale and offset if required. 22nd Nov 2010

% Ryan Eastman
% November 2010
% report issues to: rmeast@atmos.washington.edu

% tested with NCAR/NCEP reanalysis files only - so far

fname = find_file(fname);

% open the file
ncid = netcdf.open(fname,'nowrite');

% estract information on:

% numdims - # of dimensions
% numvars - # of variables
% numglobalatts - # number of global attributes
% unlimdimid - ID of dimension defined with unlimited length
[numdims numvars numglobalatts unlimdimid] = netcdf.inq(ncid);

% dimensions - struct form dims.[Dimension Name] = [Dimension Length]

for i = 0:numdims-1,
    
    % return dimension name and length
    [dname dlength] = netcdf.inqDim(ncid,i);
    
    % assign struct: dims.[Dimension Name] = [Dimension length]
    eval(['dims.',dname,' = dlength;']);
    
end;

% attributes - per variable, put in struct form: varatts.[Variable
% Name].[Attribute Name] = [Attribute Value]

for i = 0:numvars-1,
    
    % return variable name, type, ID's of dimensions and # of attributes
    [vname xtype dimids natts] = netcdf.inqVar(ncid,i);
    
    % regurn variable ID
    [varid] = netcdf.inqVarID(ncid,vname);
    
    % cycle through # of attributes per variable
    for j = 0:natts-1,
        
        % attribute name
        [atname] = netcdf.inqAttName(ncid,varid,j);

        % attribute value
        at = netcdf.getAtt(ncid,varid,atname);
        
        % ugly fix for attributes starting with '_' that break the eval line
        % specifically _FillValue and FillValue_
        if atname(1) == '_', atname = atname(2:end); end;
        if atname(end) == '_', atname = atname(1:end-1); end;
        
        % assign struct: varatts.[Variable Name].[Attribute Name] =
        % [Attribute Value]
        eval(['varatts.',vname,'.',atname,' = at;']);

    end;
    
end;

% some redundancy here for simplicity...

% variable values in struct form: vars.[Variable Name] = [Variable Value]

for i = 0:numvars-1,
    
    % return variable name and length
    [vname vlength] = netcdf.inqVar(ncid,i);
    
    % return variable ID
    vid = netcdf.inqVarID(ncid,vname);
    
    % return variable value
    vval = netcdf.getVar(ncid,vid);
    
    % Replace missing values with NaN
    s = 'missing_value';   
    eval(['ind = isfield(varatts.', vname,', s);']);
    
    if ind
        eval(['vval = fillnan(double(vval),varatts.',vname,'.missing_value);']);
    end
    
    % Replace filled values with NaN
    s = 'FillValue';   
    eval(['ind = isfield(varatts.', vname,', s);']);
    
    if ind
        eval(['if isnumeric(varatts.',vname,'.FillValue), vval = fillnan(double(vval),varatts.',vname,'.FillValue);'...
              'else vval = fillnan(double(vval),str2num(varatts.',vname,'.FillValue)); end']);
    end
    
    % Scale if applicable
    s = 'scale_factor';
    eval(['ind = isfield(varatts.', vname,', s);']);
    if ind
       eval(['vval = vval*varatts.',vname,'.scale_factor;']);
    end
    
    % Offset if applicable
    s = 'add_offset';
    eval(['ind = isfield(varatts.', vname,', s);']);
    if ind
       eval(['vval = vval + varatts.',vname,'.add_offset;']);
    end
    
    % assign struct: vars.[Variable Name] = [Variable Value]
    eval(['vars.',vname,' = vval;']);
end;

% close netcdf file
netcdf.close(ncid);