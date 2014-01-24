function [out] = readcolorsheme()
    fid = 'colorsheme18v3.txt';
    fopen(fid,'r');
    out = textread(fid);
    %colorsheme = timeseries(temp{1,1},temp{1,2},temp{1,3});
    fclose('all');