function [INV] = inventory(indir0,td_to_get)
% 5/21/2013  Parker MacCready
%
% this makes an inventory of times in processed NCOM (or other)
% files.  This will then returns a list of filenames and indices
% corresponding to a desired time record.
%
% NOTE: NCOM files are stored with a vector "time" that is time in seconds
% from the start of the year, at the start of each day, and the vectors are
% of length 365 or 366.  So for a non-leap year the first time is 0, and
% the last is 364*86400 (the start of day 365).
%
%
% INPUT:
% * indir0 is the place where all the NCOM processed folders (e.g. "pro2006")
% are stored
% * td_to_get is a vector of datenum times for which INV will give the
% information needed for make_clim.m to find them
%
% OUTPUT: each row in the one of the output variables is for a requested
% time (as specified by td_to_get) so INV.td should equal td_to_get, except
% perhaps for some rounding done by the dsearchn call
%
% * INV.td = datenum time
% * INV.indir_vec = directory name (only different each year)
% * INV.index_vec = time indices in a given year-long file

yvec = 1980:2020; % the list of years to look for
tag = 'pro'; % tag used in the "processed" directory names (like pro2006/)
D = dir(indir0); % get s structure of the directory contents

% first make a list of the names of existing data directories
kk = 0; % a counter
for ii = 1:length(yvec)
    y = yvec(ii);
    ys = num2str(y);
    name = [tag,ys];
    for jj = 1:length(D);
        iy = strcmp(name,D(jj).name);
        if iy
            kk = kk+1;
            DD.name{kk} = name;
            DD.y{kk} = y;
        end
    end
end

% now look in each directory and find the time vector of the data, using
% ssh.nc as representative of all
for ii = 1:length(DD.name)
    indir = DD.name{ii};
    infile = 'ssh.nc';
    fn = [indir0,indir,'/',infile];
    t = nc_varget(fn,'time'); % time in seconds since the start of the year
    td = datenum(DD.y{ii},1,1,0,0,0) + t/86400; % and a datenum time
    DD.t{ii} = t;
    DD.td{ii} = td;
end

% now string these out to make searchable vectors
% initialize vectors (they will be row vectors)
td_vec = []; % datenum time
indir_vec = []; % directory name
index_vec = []; % time indices in a given year-long file
for ii = 1:length(DD.name)
    td = DD.td{ii}; % will be a vector
    td_vec = [td_vec; td];
    indir_vec = [indir_vec; repmat(DD.name{ii},length(td),1)];
    index_vec = [index_vec; [1:length(td)]'];
end

% check for out of range values
if td_to_get(1)<td_vec(1) | td_to_get(end)>td_vec(end)
    disp('*** requested times out of range !!! ***')
    return
end

% indices into the searchable vectors for the requested times
isv = dsearchn(td_vec,td_to_get(:));

% store the results
INV.td = td_vec(isv);
INV.indir = indir_vec(isv,:);
INV.index = index_vec(isv);

