function theResult = ncload2(theNetCDFFile, varargin)

% NCLOAD2: as ncload, but missing_value and _FillValue are replaced by NAN,
% variable ranges can be specified, and data is multiplied by scale_factor 
% if present.
% Example: ncload2('test.nc','T(:,:,12:21,24:73)','S')
% Disadvantages: Takes about twice as long to run as ncload because
% variables must actually be evaluated to check values and replace some
% with NANs; could probably be sped up if coded better. It can still be
% quite a bit faster than ncload if you're only reading in partial ranges
% of large variables.
% Ian Eisenman, 2006
%
% ncload -- Load NetCDF variables.
%  ncload('theNetCDFFile', 'var1', 'var2', ...) loads the
%   given variables of 'theNetCDFFile' into the Matlab
%   workspace of the "caller" of this routine.  If no names
%   are given, all variables are loaded.  The names of the
%   loaded variables are returned or assigned to "ans".
%   No attributes are loaded.
 
% Copyright (C) 1997 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 18-Aug-1997 10:13:57.

if nargin < 1, help(mfilename), return, end

result = [];
if nargout > 0, theResult = result; end

f = netcdf(theNetCDFFile, 'nowrite');
if isempty(f), return, end

if isempty(varargin), varargin = ncnames(var(f)); end

for i = 1:length(varargin)
   if ~isstr(varargin{i}), varargin{i} = inputname(i+1); end
   % check if range is specified. if it is, strip it from varargin{i} and
   % load variable only in range.
   if sum(varargin{i}=='(') > 0
       g=find(varargin{i}=='(');
       var_range=varargin{i}(g(1):end);
       varargin{i}=varargin{i}(1:g(1)-1);
       var0=eval(['f{varargin{i}}' var_range]);
   else
       var0=f{varargin{i}}(:);
   end
   % replace missing_value and _FillValue with NAN
   mv=f{varargin{i}}.FillValue_(:); if ~isempty(mv), var0(var0==mv)=nan; end
   mv=f{varargin{i}}.missing_value(:); if ~isempty(mv), var0(var0==mv)=nan; end
   % multiply by scale_factor
   sf=f{varargin{i}}.scale_factor(:); if ~isempty(sf), var0=var0*sf; end
   assignin('caller', varargin{i}, var0)
end

result = varargin;

close(f)

if nargout > 0
   theResult = result
else
   ncans(result)
end
