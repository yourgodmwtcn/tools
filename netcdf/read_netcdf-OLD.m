function data = read_netcdf(filename,varargin)
%data=READ_NETCDF(filename[,'flipdim']) Extracts netcdf file to structure
%   This function extracts a netcdf file into a structure which is
%   returned to the command line.  You need only specify a filename of a
%   netCDF file.  This version uses the java class to read data directly
%   from the netCDF file.
%
%   Options:
%   READ_NETCDF(filename,'flipdim') This will automatically permute the
%   arrays read from the file in reverse order.
%
%   READ_NETCDF(filename,'strip') This will strip non-alphanumerica
%   characters from element names to allow them to be used as structure
%   elements.
%
%   Note:
%   Some versions of MATLAB will require a local copy of the netCDF .jar
%   files.  To temporarily add these use:
%   javaaddpath('ftp://ftp.unidata.ucar.edu/pub/netcdf-java/v4.1/netcdf-4.1.jar');
%   Or you can download the files locally and add them to your classpath.txt
%   file.  You can find your classpath.txt file using the which command.
%   which classpath.txt
%
%   Example
%       data=read_netcdf('input.nc');
%
%   Version 1.1
%   Maintained by: Samuel Lazerson (lazerson@pppl.gov)
%   Date  10/07/2010

% Allow the user to pass some variables.
flipdim=0;
strip=0;
if nargin>1
    for i=2:nargin
        switch varargin{i-1}
            case 'flipdim'
                flipdim=1;
            case 'strip'
                strip=1;
        end
    end
end
% Add the path to your downloaded JAR file and restart matlab.
import ucar.nc2.*;
%%%%Open the File%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
netcdfdata=NetcdfFile.open(filename);
% Get information on number of elements
ndimen=netcdfdata.getDimensions.size;
nvars=netcdfdata.getVariables.size;
ngatts=netcdfdata.getGlobalAttributes.size;
unlimidimid=netcdfdata.getUnlimitedDimension;
% Get global attributes
if ngatts >0
    for i=0:ngatts-1
        att=netcdfdata.getGlobalAttributes.get(i);
        attname=char(att.getName);
        if strip, attname(~isstrprop(attname,'alphanum'))=''; end;
        % We handle multiple data types
        if (att.getDataType.isString)
            data.(attname)=char(att.getStringValue);
        elseif (att.getDataType.isNumeric)
            data.(attname)=att.getValue;
        end
    end
end
% Get dimensions
if ndimen >0
    for i=0:ndimen-1
        dim=netcdfdata.getDimensions.get(i);
        dimname=char(dim.getName);
        if strip, dimname(~isstrprop(dimname,'alphanum'))=''; end;
        data.(dimname)=dim.getLength;
    end
end
% Get Variables
if nvars >0
    for i=0:nvars-1
        var=netcdfdata.getVariables.get(i);
        varname=char(var.getName);
        if strip, varname(~isstrprop(varname,'alphanum'))=''; end;
        natts=var.getAttributes.size;
        if (var.getDataType.isString)
            try
                data.(varname)=char(var.getStringValue);
            catch exception
                try
                    data.(varname)=char(var.readScalarString);
                catch excpetion
                    disp(['-----Error reading: ' varname]);
                end
            end
        elseif (var.getDataType.isNumeric)
            if var.getSize == 1
                data.(varname)=var.read.copyTo1DJavaArray;
            else
                if flipdim
                    temp=var.read.copyToNDJavaArray;
                    temp=permute(temp,ndims(temp):-1:1);
                    data.(varname)=temp;
                else
                    data.(varname)=var.read.copyToNDJavaArray;
                end
            end
        end
        if natts > 0
            for j=0:natts-1
                att=var.getAttributes.get(j);
                attname=char(att.getName);
                if strip, attname(~isstrprop(attname,'alphanum'))=''; end;
                varattname=strcat(varname,'_',attname);
                if (att.getDataType.isString)
                    data.(varattname)=char(att.getStringValue);
                elseif (att.getDataType.isNumeric)
                    data.(varattname)=att.getValue;
                end
            end
        end
    end
end
% Close the file
netcdfdata.close;
end

