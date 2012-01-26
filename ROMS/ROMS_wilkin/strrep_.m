function s = strrep_(s)
% $Id: strrep_.m 358 2008-04-07 14:15:03Z zhang $
% Replace all underscore characters in a string with \_ so that the
% underscore is not interpretted as a TeX subscript instruction.
% This function is used by several roms_*view routines so that ROMS
% filenames are not corrupted in the title string
%
% John Wilkin
s = strrep(s,'_','\_');
