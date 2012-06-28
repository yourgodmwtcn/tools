function roms_averageFileSet(seriesIn, ncn, filename, weights, outdir)

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

% loop through variables inside
info = nc_info(tmp);
vars = info.Dataset;
for vi=1:length(vars)
	if length(vars(vi).Size) >= 2, disp(['    ' vars(vi).Name]); end
	if vars(vi).Nctype ~= nc_char
		% loop through files and make a weighted average of that variable
		A = weights(1) * nc_varget(inNames{1}, vars(vi).Name);
        for ni = 2:length(ncn)
            A = A + (weights(ni) * nc_varget(inNames{ni}, vars(vi).Name));
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
		nc_varput(tmp, vars(vi).Name, A);
	end
end

% if the output filename is a relative path, prepend the series' dirname
% commented out by PM 3/24/2011
% if filename(1) ~= '/' & filename(1) ~= '~'
% 	filename = [seriesIn.dirname filename];
% end

% name the output file correctly to indicate that we're done
system(['mv ' tmp ' ' filename]);
dt = toc;
disp(['    done: created ' filename,' in ',num2str(dt),' sec']);
