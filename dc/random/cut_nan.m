% snips out non-nan sections and concatenates them
%   [out] = cut_nan(in)

function [out] = cut_nan(in)
    
    out = in(~isnan(in));