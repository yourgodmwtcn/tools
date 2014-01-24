function [data,coords,grid] = roms_extractFromSeries(series, varname, t, theType, varargin);

% [data,coords] = roms_extractFromSeries(series, varname, t, ...
% [data,coords,grid] = ...
%
% general routine for extracting data from a ROMS netcdf file series
% (e.g., ocean_his_*.nc) in data units. _series_ is a structure created by 
% roms_createSeriesDef.m. There should be no reason for the user to call
% this directly: use roms_extract.m. See roms_extract.m for full syntax.
%
% neil banas, 2009-2011


if ~isstruct(series) | ~isfield(series,'ncn') | ~isfield(series,'nctime')
	error('need to supply a seriesDef with ncn and nctime defined.');
elseif isempty(series.ncn) | isempty(series.nctime)
	error('file series timebase is empty.');
end

if length(t) > 1 % if t is a vector, recurse over each time in the vector

	outputsInitialized = false;
	for n = 1:length(t)
		% extract one time slice
		disp(['    extracting ' varname ' ' num2str(n) ' / ' num2str(length(t))]);
		if n==1
			[data1, coords1, grid] = roms_extractFromSeries(series, varname, t(n), theType, varargin{:});
		else
			[data1, coords1] = roms_extractFromSeries(series, varname, t(n), theType, varargin{:}, 'grid', grid);
		end
		if ~isempty(data1) 
			if ~outputsInitialized
				outputsInitialized = true;
				data = repmat(nan,[length(t) size(data1)]);
				coords.tm = data;
				if isfield(coords1,'zm')
					coords.zm = data;
				end
				coords.ym = data;
				coords.xm = data;
				L = length(data1(:));
			end
			% place the time slice in the output variables
			data(n,:) = reshape(data1, [1 L]);
			coords.tm(n,:) = repmat(t(n), [1 L]);
			if isfield(coords,'zm')
				coords.zm(n,:) = reshape(coords1.zm, [1 L]);
			end
			coords.ym(n,:) = reshape(coords1.ym, [1 L]);
			coords.xm(n,:) = reshape(coords1.xm, [1 L]);
		end
	end
		
else % t is scalar

	% file numbers bracketing t
	n = interp1(series.nctime, series.ncn, t);
	if isnan(n)
		error(['can''t find file numbers to go with ' datestr(t)]);
	else		
		n0 = max(floor(n),series.ncn(1));
		file0 = roms_filename([series.dirname series.basename],n0);
		n1 = min(ceil(n),series.ncn(end));
		file1 = roms_filename([series.dirname series.basename],n1);
		% do the extraction
		[data,coords,grid] = roms_extractFromFile(file0, varname, theType, varargin{:});
		% interpolate between frames if necessary
		if n0 ~= n1
			[data1,coords1] = roms_extractFromFile(file1, varname, theType, varargin{:});
			fr = (n-n0)/(n1-n0);
			data = data + (data1-data).*fr;
			if isfield(coords,'zm')
				coords.zm = coords.zm + (coords1.zm-coords.zm).*fr;
			end
		end
	end
	
end

