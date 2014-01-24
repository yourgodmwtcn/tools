function [data,coords] = roms_extractFromVariable(data0, coords0, theType, varargin);

% [data,coords] = ...
%     roms_extractFromFullVariable(full3Dvardata, fullcoords, 'point', y, x);
%     roms_extractFromFullVariable(full4Dvardata, fullcoords, 'zslice', z);
%                                                        ..., 'depthslice', depth);
%                                                        ..., 'depthaverage', [mindepth maxdepth]);
%                                                        ..., 'depthintegral', [mindepth maxdepth]);
%                                                        ..., 'profile', y, x);
%                                                        ..., 'point', z, y, x);
%
% extracts various shapes from full variables read by roms_extract. Useful for manually
% chaining together a bunch of extractions.
%
% works on single-frame extractions only: "full3Dvardata" refers to the case of lat-lon fields
% like zeta, and "full4Dvardata" refers to the case of depth-lat-lon fields like salt.
%
% neil banas, 2013


if ndims(data0)==3
% 3D variables (time-lat-lon, one time frame only) ---------------------------------------
	x2 = coords0.xm;
	y2 = coords0.ym;
	switch theType
		case {'profile','profiles','point','points'} % ----- 3D, points
			ym = varargin{1}(:);
			xm = varargin{2}(:);
			L = max([length(ym) length(xm)]);
			if length(xm)==1, xm = repmat(xm,[L 1]); end
			if length(ym)==1, ym = repmat(ym,[L 1]); end
			data = interp2(x2,y2,data0,xm,ym);
			
		otherwise
			error(['don''t know how to extract ''' theType ''' for 3D variables.']);
	end
	
else
% 4D variables (time-depth-lat-lon, one time frame only) ---------------------------------


% *** rewrite starting here
		case {'zslice','zslices','zSlice','zSlices'} % ----- 4D, slice at one or more z values (from MSL)
			z = varargin{1};
			zm = repmat(z(:),[1 J I]);
			ym = repmat(reshape(y2,[1 J I]),[length(z(:)) 1 1]);
			xm = repmat(reshape(x2,[1 J I]),[length(z(:)) 1 1]);
			zetam = repmat(reshape(zeta2,[1 J I]),[length(z(:)) 1 1]);
			Hm = repmat(reshape(H2,[1 J I]),[length(z(:)) 1 1]);
			csm = (zm - zetam) ./ (zetam + Hm);
			bad = isnan(csm) | csm > 0 | csm < -1;
			csm(isnan(csm)) = 0;
			data3 = squeeze(nc_varget(filename,varname));
			x3 = repmat(reshape(x2, [1 J I]), [K 1 1]);
			y3 = repmat(reshape(y2, [1 J I]), [K 1 1]);
			cs3 = repmat(cs(:),[1 J I]);
			data = interpn(cs3,y3,x3,data3,csm,ym,xm);
			data(bad) = nan;
			
		case {'depthslice','depthslices','depthSlice','depthSlices'} % ----- 4D, slice at one or more depth values (from surface)
			depth = varargin{1};
			zetam = repmat(reshape(zeta2,[1 J I]),[length(depth(:)) 1 1]);
			zm = repmat(depth(:),[1 J I]) + zetam; % this is the only difference between depthslice and zslice
			ym = repmat(reshape(y2,[1 J I]),[length(depth(:)) 1 1]);
			xm = repmat(reshape(x2,[1 J I]),[length(depth(:)) 1 1]);
			Hm = repmat(reshape(H2,[1 J I]),[length(depth(:)) 1 1]);
			csm = (zm - zetam) ./ (zetam + Hm);
			bad = isnan(csm) | csm > 0 | csm < -1;
			csm(isnan(csm)) = 0;
			data3 = squeeze(nc_varget(filename,varname));
			x3 = repmat(reshape(x2, [1 J I]), [K 1 1]);
			y3 = repmat(reshape(y2, [1 J I]), [K 1 1]);
			cs3 = repmat(cs(:),[1 J I]);
			data = interpn(cs3,y3,x3,data3,csm,ym,xm);
			data(bad) = nan;
		
		case {'depthaverage','depthAverage'} % ----- 4D, average over some depth range
			ym = y2; % this can't be vectorized; output is always 2D
			xm = x2;
			data3 = squeeze(nc_varget(filename,varname));
			[foo,data] = roms_depthIntegrate(data3, G.cs, G.csw, H2, zeta2, varargin{1});
			
		case {'depthintegral','depthIntegral'} % ----- 4D, integral over some depth range
			ym = y2; % this can't be vectorized; output is always 2D
			xm = x2;
			data3 = squeeze(nc_varget(filename,varname));
			[data,foo] = roms_depthIntegrate(data3, G.cs, G.csw, H2, zeta2, varargin{1});

		case {'profile','profiles'} % ----- 3D, vertical profiles at one or more (y,x) pairs
			y = varargin{1}(:);
			x = varargin{2}(:);
			if length(x)==1 & length(y)>1
				x = repmat(x,size(y));
			elseif length(x)>1 & length(y)==1
				y = repmat(y,size(x));
			end
			L = length(x);
			cs_extrap = unique([-1; sort(cs); 0]); % extend to top & bottom of water column
			Km = length(cs_extrap);
			ym = repmat(reshape(y,[1 L]),[Km 1]);
			xm = repmat(reshape(x,[1 L]),[Km 1]);
			zetam = interp2(x2,y2,zeta2,xm,ym);
			Hm = interp2(x2,y2,H2,xm,ym);
			csm = repmat(cs_extrap,[1 L]);
			zm = zetam + csm .* (zetam + Hm);
			data3 = squeeze(nc_varget(filename,varname));
			x3 = repmat(reshape(x2, [1 J I]), [K 1 1]);
			y3 = repmat(reshape(y2, [1 J I]), [K 1 1]);
			cs3 = repmat(cs(:),[1 J I]);
			data = interpn(cs3,y3,x3,data3,csm,ym,xm);
		
		case {'point','points'} % ----- 4D, data at one or more (z,y,x) triplets
			zm = varargin{1}(:);
			ym = varargin{2}(:);
			xm = varargin{3}(:);
			L = max([length(xm) length(ym) length(zm)]);
			if length(xm)==1, xm = repmat(xm,[L 1]); end
			if length(ym)==1, ym = repmat(ym,[L 1]); end
			if length(zm)==1, zm = repmat(zm,[L 1]); end
			zetam = interp2(x2,y2,zeta2,xm,ym);
			Hm = interp2(x2,y2,H2,xm,ym);
			csm = (zm - zetam) ./ (zetam + Hm);
			bad = isnan(csm) | csm > 0 | csm < -1;
			csm(isnan(csm)) = 0;
			data3 = squeeze(nc_varget(filename,varname));
			x3 = repmat(reshape(x2, [1 J I]), [K 1 1]);
			y3 = repmat(reshape(y2, [1 J I]), [K 1 1]);
			cs3 = repmat(cs(:),[1 J I]);
			data = interpn(cs3,y3,x3,data3,csm,ym,xm);
			data(bad) = nan;
			
		otherwise
			error(['don''t know how to extract ''' theType ''' for 4D variables.']);
	end
end