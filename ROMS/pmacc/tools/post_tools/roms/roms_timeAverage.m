function Sout = roms_timeAverage(Sin, period, startTime, nAvgs, suffix);

% seriesDefOut = roms_timeAverage(seriesDefIn, period, startTime, nAvgs, suffix);
% seriesDefOut = roms_timeAverage(seriesDefIn, period, suffix);
%
% chunk-averages a series of .nc files into a new series.
%
% default startTime is the start of the time series; default # of averages is Inf
% (i.e., all of them).
%
% neil banas feb 2009

if nargin<4
	suffix = startFrame;
	startFrame = 1;
	nAvgs = Inf;
end

if length(Sin.ncn)<2, error('no file series specified.'); end
dt = Sin.nctime(2) - Sin.nctime(1);
nPerAvg = round(period / dt);

iAvg = 1;
n0 = round(interp1(Sin.nctime,Sin.ncn,startTime);
nn = n0 : n0+nPerAvg-1;
while nn(end) <= Sin.ncn(end)  &  iAvg <= nAvgs
	outname = roms_filename([Sin.dirname Sin.basename suffix '_'], iAvg);
	roms_averageFileSet(Sin, nn, outname);
	iAvg = iAvg+1;
	n0 = n0 + nPerAvg;
	nn = n0 : n0+nPerAvg-1;
end

Sout = roms_createSeriesDef(Sin.dirname, [Sin.basename suffix '_']);