function xf = llpa(x,nf)

% LLPA low pass filter x into xf
% Call as 
%	xf = llpa(x,nf)
% Data come out at 6 hour intervals
% nf = first central point of data from input series
%	(nf >=60)
% Input data are hourly

% On #2, the mfile reads the filter weights (the. DAT file). 
% You give it hourly data and it automatically returns filtered data at 6 hour intervals.
% I like to run the time series of time (as decimal days) through as well so you know what times 
% correspond to what sea levels (or whatever). You will probably want to skip a few data points 
% at the beginning so the times come out even. With the 6 hourly data you can either daily average those, 
% or simply just sample one point per day. The HPP (half power point) is something like 46 hours, 
% so you do not really lose information if you do that.

load LPW.DAT
lpw = LPW;
[nff,m] = size(lpw);
ff = zeros(2*nff-1,1);
ff(1:nff) = lpw(nff:-1:1);
ff(nff:((2*nff)-1))= lpw;

iff = 1;
nl = nf - 59;
nh = nf + 59;
[n,m] = size(x);

while nh <= n
	y = x(nl:nh,:);
	for j = 1:m
		xf(iff,j) = ff'*y(:,j);
	end
	iff = iff + 1;
	nf = nf + 6;
	nl = nf - 59;
	nh = nf + 59;
end

