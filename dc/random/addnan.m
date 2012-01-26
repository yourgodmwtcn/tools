% addnan(var,val) replaces all values in (var) > val with NaN

function [A] = addnan(var,val)
    aa = (var) > val;
    A = var;
    A(aa) = NaN;