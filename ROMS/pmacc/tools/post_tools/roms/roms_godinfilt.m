function Sout = roms_godinfilt(Sin, outputTimebase, suffix,outdir,istart)

% seriesDefOut = roms_godinfilt(seriesDefIn, outputTimebase, suffix, outdir, istart);
% seriesDefOut = roms_godinfilt(seriesDefIn, suffix);
%
% filters all the variables in a file series using a 24-24-25 hr window and saves
% a new series of the results (Godin filter).
%
% NOTE: Sin must be hourly saves and default is the 24-24-25 hr window, with one output file
% per day, centered on noon.
%
% neil banas feb 2009, modified by DAS, Oct. 2010
% edited by SNG 4/16/2012 to include istart as an option in order to start the code on an output file other than 1
% and included outdir in the file structure
% SNG 4/16/2012 updated t0 and t1 to properly represent the symmetrical
% filter and the weights calculation to account for roundoff errors in
% Sin.nctime
% SNG 4/17/2012 checked the output from this code with a low pass godin
% filter on ssh at an individual mooring and they check out perfectly
%------------------------------------------------------------------------

if nargin<3
	suffix = outputTimebase;
	outputTimebase = [];
    istart = 1;
end

% Godin filter is 71 hours long
filterWindow = 71/24; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build weights for Godin filter- start with 24 hr. filter
 filter24 = ones(24,1); 
 filter24 = filter24 ./ sum(filter24);
 
% then build 25 hr filter
 filter25 = ones(25,1); 
 filter25 = filter25 ./ sum(filter25);
 
% covolve filters together, works because conv is associative
 temp_filter = conv(filter24, filter24);
 filter = conv(temp_filter, filter25);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    %starting time = central time - half of the window length
    %the -1/24 is included because this is a symmetrical filter 
	t0 = ti - (filterWindow-1/24)/2;
	t1 = ti + (filterWindow-1/24)/2;
	f = find(Sin.nctime >= ti-filterWindow/2 & Sin.nctime <= ti+filterWindow/2);
    %if the time stamps of the files exactly match the desired time stamps
    %(i.e. you have hourly files) then choose weights = filter
    %these statements are needed to correct round-off errors in Sin.nctime
    if length(f) == length(filter)
        if linspace(t0,t1,length(filter))-Sin.nctime(f)<10^-9
            weights = filter';
        else
            weights = interp1(linspace(t0,t1,length(filter)), filter, Sin.nctime(f));
        end
    else
        weights = interp1(linspace(t0,t1,length(filter)), filter, Sin.nctime(f));
    end
    nn = Sin.ncn(f);
    outname = roms_filename([outdir Sin.basename suffix '_'], iAvg); % SNG to match PM roms_hanning.m edit
    roms_averageFileSet(Sin, nn, outname, weights, outdir);
    iAvg = iAvg+1;
end

Sout = roms_createSeriesDef(Sin.dirname, [Sin.basename suffix '_']);	