function[yf]=yearfrac(num)
%YEARFRAC   Convert date from 'datenum' format to 'year.fraction'.
%  
%   YF=YEARFRAC(NUM) where NUM is an array of dates in Matlab's 'datenum'
%   format, returns the fraction of the year at each date.
%
%   The actual number of days in each year is used, that is, including 
%   leap years.
%
%   YR=YEARFRAC(NUM) where NUM is a cell array of numeric arrays, also
%   works.  YR is thena cell array of the same size as NUM.
%  
%   See also YF2NUM, DATENUM, DATEVEC
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 1998--2012 J.M. Lilly --- type 'help jlab_license' for details        
  
if ~iscell(num)
    if strcmp(num, '--t')
        yf2num('--t'),return
    end  
end

if iscell(num)
    for i=1:length(num)
        yf{i,1}=yearfrac1(num{i});
    end
else
    yf=yearfrac1(num);
end

function[yf]=yearfrac1(num)
yf=[];
if ~isempty(num)
    index=find(isnan(num));
    if ~isempty(index)
        num(index)=0;
    end

    [y,mo,d,h,mi,s] = datevec(num);

    %Number of days in this year?
    nd=datenum(y,12,31)-datenum(y-1,12,31);


    na=datenum(y,mo,d,h,mi,s)-datenum(y-1,12,31)-1;
    %The minus one is because for Jan 1, I add nothing to year.fraction

    yf=y+na./nd;

    if ~isempty(index)
      yf(index)=nan;
    end
end

  
  
