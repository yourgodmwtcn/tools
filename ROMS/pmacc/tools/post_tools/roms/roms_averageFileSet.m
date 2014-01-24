function roms_averageFileSet_native(seriesIn, ncn, filename, weights, outdir)

% roms_averageFileSet(seriesIn, frameNums, outputFilename);
%                                                    ...,weights);
%															..., outdir);
%
% averages a set of netcdf files into a single netcdf file with the same format.
% looks for files with numbers _ncn_, which should be a subset of seriesIn.ncn.
% All variables are averaged, except for those of type nc_char.
% 
% a vector of weights the same length as frameNums can be passed, in order to
% use this routine for filtering.
%
% if outdir is not given, uses seriesIn.dirname: i.e., by default makes the
% new files in the same place the input files are found, but doesn't have to.
%
% neil banas feb 2009

% with some editing by PM 3/24/2011 to reshape arrays so they match those
% in the original history file

% updated 10/6/2011 SNG to include biological variables NO3, phytoplankton,
% zooplankton, detritus, and oxygen
% updated 1/24/2012 SNG to include dye variables
% SNG 18 April 2012 incorporated native matlab netcdf calls rather than
% snctools. The speed up was only about 5% though so not much. Otherwise,
% this code is identical to roms_averageFileSet.m and I checked the output
% results and they are identical

if nargin < 4
	weights = ones(size(frameNums)) ./ length(frameNums);
end

if nargin < 5
	outdir = seriesIn.dirname;
end

% assemble filenames
inNames = cell(length(ncn),1);
for ni = 1:length(ncn)
	n = ncn(ni);
	inNames{ni} = roms_filename([seriesIn.dirname seriesIn.basename], n);
	if ~exist(inNames{ni},'file'), error(['can''t find ' inNames{ni}]); end
end

tic
disp(' ');
disp(['averaging ' seriesIn.dirname seriesIn.basename ' #' num2str(ncn(1)) ' to ' num2str(ncn(end))]);

% duplicate the first input file to create a template for the output
%tmp = [seriesIn.dirname 'temp.nc'];
%SNG temporarily changed to last file because some of the salish_2005_1
%output files have Akv and Aks and some do not!
tmp = [outdir 'temp.nc'];
system(['cp ' inNames{end} ' ' tmp]);

% %I tested the below, opening all files first and there was no difference in
% %speed! SNG 19 April 2012
% %open all netcdf files to be filtered
% ncid = zeros(1,lenght(ncn));
% for ni = 1:length(ncn)
%     ncid(ni) = netcdf.open(inNames{ni},'NC_NOWRITE');
% end

% loop through variables inside
info = nc_info(tmp);
vars = info.Dataset;
for vi=1:length(vars)
	if length(vars(vi).Size) >= 2, disp(['    ' vars(vi).Name]); end
	if vars(vi).Nctype ~= nc_char
		% loop through files and make a weighted average of that variable
		ncid = netcdf.open(inNames{1},'NC_NOWRITE');
        var_id = netcdf.inqVarID(ncid,vars(vi).Name);
        A = weights(1) * netcdf.getVar(ncid,var_id,'double');
        netcdf.close(ncid);
        for ni = 2:length(ncn)
            ncid = netcdf.open(inNames{ni},'NC_NOWRITE');
            var_id = netcdf.inqVarID(ncid,vars(vi).Name);
            A = A + (weights(ni) * netcdf.getVar(ncid,var_id,'double'));
            netcdf.close(ncid);
            % if all files are opened ahead of time, use below 2 lines
            % rather than above 4, this is related to my time test
            %var_id = netcdf.inqVarID(ncid(ni),vars(vi).Name);
            %A = A + (weights(ni) * netcdf.getVar(ncid(ni),var_id,'double'));
        end
        % PM 3/24/2011 I had to add the reshape statement below in order to
        % get the nc_varput call to work on my macbook.  Without it the
        % size of temp.nc grew to many, many GB!
        if sum(strcmp(vars(vi).Name,{'zeta';'ubar';'vbar';'u';'v';'w'; ...
                'temp';'salt';'AKv';'AKs';'shflux';'ssflux';'sustr';'svstr'; ...
                'NO3';'phytoplankton';'zooplankton';'detritus';'oxygen';'dye_'})) ...
            || sum(strncmp(vars(vi).Name,'dye_',4))
            A = reshape(A,[1 size(A)]); % pm 3/24/2011
            %disp([vars(vi).Name,': size of A = ',num2str(size(A))]);
        end
        % write the average into the output file
        ncid = netcdf.open(tmp,'NC_WRITE');
        var_id = netcdf.inqVarID(ncid,vars(vi).Name);
        netcdf.putVar(ncid,var_id,A);
        netcdf.close(ncid);
	end
end

% %also related to my test, if they were all opened at once, they all need to
% %be closed
% %close all netcdf files
% for ni = 1:length(ncn)
%     netcdf.close(ncid(ni));
% end

% if the output filename is a relative path, prepend the series' dirname
% commented out by PM 3/24/2011
% if filename(1) ~= '/' & filename(1) ~= '~'
% 	filename = [seriesIn.dirname filename];
% end

% name the output file correctly to indicate that we're done
system(['mv ' tmp ' ' filename]);
dt = toc;
disp(['    done: created ' filename,' in ',num2str(dt),' sec']);
