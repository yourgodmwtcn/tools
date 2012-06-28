function Sout = roms_hanning(Sin, filterWindow, outputTimebase, suffix, outdir,istart)

% seriesDefOut = roms_hanning(seriesDefIn, filterWindow, outputTimebase, suffix, outdir);
% seriesDefOut = roms_hanning(seriesDefIn, suffix, outdir);
%
% filters all the variables in a file series using a hanning window and saves
% a new series of the results.
%
% default is a 40 hr window (i.e., filterWindow = 40/24) and one output file
% per day, centered on noon.
%
% neil banas feb 2009

% edited by PM 3/24/2011
% edited by SNG 6/20/2011 to include istart as an option in order to start the code on an output file other than 1

if nargin<5
	suffix = filterWindow;
	filterWindow = [];
	outputTimebase = [];	
	istart = 1;
end
if isempty(filterWindow)
	filterWindow = 40/24;
end

if length(Sin.ncn)<2, error('no file series specified.'); end
t0 = Sin.nctime(1) + filterWindow/2;
t1 = Sin.nctime(end) - filterWindow/2;
if isempty(outputTimebase)
	% default: one file per day at noon
	outputTimebase = ceil(t0)+0.5 : t1;
else
	% make sure the output times requested fall within the file series,
	% with a margin for half the filter window on each end
	outputTimebase = outputTimebase(outputTimebase >= t0 & outputTimebase <= t1);
end

iAvg = istart; %SNG changed this line and the following to start at iAvg = start so you can run this code starting at an output file other than 1
for i = iAvg:length(outputTimebase)
	ti = outputTimebase(i);
	t0 = ti - filterWindow/2;
	t1 = ti + filterWindow/2;
	f = find(Sin.nctime >= ti-filterWindow/2 & Sin.nctime <= ti+filterWindow/2);
	weights = interp1(linspace(t0,t1), 1+cos(linspace(-pi,pi)), Sin.nctime(f));
	weights = weights./sum(weights);
	nn = Sin.ncn(f);
	outname = roms_filename([outdir Sin.basename suffix '_'], iAvg); % PM edit
	roms_averageFileSet(Sin, nn, outname, weights,outdir);
	iAvg = iAvg+1;
end

Sout = roms_createSeriesDef(Sin.dirname, [Sin.basename suffix '_']);	