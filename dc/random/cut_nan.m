% snips out non-nan sections and concatenates them

function [out] = cut_nan(in)
    
    out = in(~isnan(in));