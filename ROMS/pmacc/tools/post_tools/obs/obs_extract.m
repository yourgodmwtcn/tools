function [data, filenames] = obs_extract(files, vars, timeRange, varargin)
%------------------------------------------------------------
% [data, filenames] = obs_extract(filename, vars, timeRange
%                                         ..., 'section',x,y,range);
%                                           ...,'section',track,range)
%                                         ..., 'point',x,y);
%                                         ..., 'polygon',x,y);
%                                         ..., 'all');
%[data, filenames] = obs_extract(directory, ...
% [data, filenames] = obs_extract(cell array of files and directories, ...
%
% looks inside a netcdf file for one or more observational variables and
% returns them in a structure _data_, along with coordinate variables.
%
% if a directory is given instead of a filename, recurses over
% files inside it and concatenates the output.
%
%Only includes data from files with single dimensions, skips mooring data
%
% vars can be a string containing one variable name, or a cell array of strings,
% or 'all'.
%
% timeRange is a 2 element vector of matlab datenums or 'all'.
%
% if 'point' is used, it looks for data within 0.1 km of that location 
%
% written by C. Bassin ,D. Sutherland and N. Banas, UW, Jan 2010
%
% edited Dec 2010 by C. Bassin:
% added filenames and data.fileid as an output.
%   filenames shows files with correct dimension size that were searched 
%   data.fileid gives an integer value that is associated to filenames cell
%   aray to show which file the data originated

%------------------------------------------------------------

debug = 0;

[netcdfpaths,files]=dirwalker(files);

x=1;
for i=1:length(netcdfpaths);

    Info_file=nc_info(netcdfpaths{i});
   
    if size(Info_file.Dimension,1)==1
        new_netcdfpaths{x}=netcdfpaths{i};
        x=x+1;
    else
        if debug
        	disp(['skipping ' netcdfpaths{i} ' due to incorrect dimension size'])
        end
    end

end

if exist('new_netcdfpaths','var')==0
    data=[];
    return
end

N=length(new_netcdfpaths);
dataC=cell(N,1);


for j=1:N;
  if debug, disp(['reading data from ' new_netcdfpaths{j}]); end

  dataC{j,1}= obs_extractFromFile(new_netcdfpaths{j}, vars, timeRange, varargin{:});
    dataC{j,1}.fileid=ones(size(dataC{j,1}.z))*j;  % puts a fileid in extracted data
  dataC{j,1}.filem=1;% create an id to tell obs merge not to get rid of fileid
end 

if N>1

    data=obs_merge(dataC,'any');
else
    data= dataC{1,1};
end
data = rmfield(data, 'filem');

% get rid of nan data in coordinate variables
data=obs_omit(data,'x == nan');
data=obs_omit(data,'y == nan');
data=obs_omit(data,'z == nan');
data=obs_omit(data,'t == nan');


data.cast=obs_identifyCasts(data);

filenames=new_netcdfpaths';

end

