function [smooth] = Z_8d(data)
% 7/10/2013  Parker MacCready, based on jfilt by Jonathan Lilly 2005
% ** use ONLY with hourly data! **
%
% This applies the Austin-Barth 8 day filter to a vector,
% or a matrix in which each COLUMN is a data record.
%
% NOTE this may be different than their definition - it just returns a
% weighted average of the data over the previous 8 days from time t,
% with the weighting decaying from 1 to 1/e at t - 8 days.  There are 192
% hours in 8 days.
% 
% It returns a dataset of the same size you started with, padded whth
% NaN's at the ends.
%
    
filter = zeros(385,1);
filter(1:193) = exp(linspace(-1,0,193));
filter=filter./sum(filter);

n=length(filter);
smooth=zeros(size(data));
a=round(n/2);
N=size(data,1);
for i=1:size(data,2)
    temp=conv(data(:,i),filter);
    smooth(:,i)=temp(a:a+N-1);
    smooth(1:n,i)=nan*ones(n,1);
    smooth(N-n+1:N,i)=nan*ones(n,1);
end




