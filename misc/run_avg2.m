%run_avg2(data,N,dim) calculates an N-day running average of hourly data by convolving with a rectangular window.  
%                       The output is of the same length as the input so that the averaged
%                       data can be used with the original time vector and compared against 
%                       other (potentially un-averaged) quantities from the same dataset.
%                       Old version could work around NaNs, but this one can't. 
%                       2003, Tom Farrar, tomf@mit.edu
% Deepak, the only thing to watch out for here is that this code was built for making N-day running averages from hourly data, so, if you want to average over 10 timesteps, the second argument should be 10./24.






%Slow old version
%run_avg(data,N) calculates an N-day running average of hourly data.  
%                       The output is of the same length as the input so that the averaged
%                       data can be used with the original time vector and compared against 
%                       other (potentially un-averaged quantities) from the same dataset.
%                       Note that NaNs are put in the first and last N/2 slots, and that NaNs
%                       are acceptable in the input.  This routine is slow, but effective.  
%                       I recommend running it once in a data processing routine rather than
%                       computing the average repeatedly.  Calls nanmean.m.
%
%                       2001, Tom Farrar, tomf@mit.edu

function [swav]=run_avg2(srad,N,dim)

%N-day average
%g=N*24;
%swav(1:g/2)=NaN*ones(g/2,1);
%for n=1+g/2:length(srad)-g/2
%    swav(n)=nanmean(srad(n-g/2:n+g/2));
%end
%swav(1+length(srad)-g/2:length(srad))=NaN*ones(g/2,1);

%swav=swav';

%New, faster way:

M=N*24;
win=ones(M,1)./M;
[m,n]=size(srad);
if dim==2
  swav=conv2(1,win,srad,'same');
elseif dim==1
  swav=conv2(win,1,srad,'same');
elseif display('run_avg2.m can only operate on 2D arrays!')
  %break
end


