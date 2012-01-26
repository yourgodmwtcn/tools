% 
% Here are some routines to do Gregorian labeling of stacked time
% series plots, including vector stick plots.
% 
% Try using DEMOTPLT.M to see how TIMEPLT and STACKLBL work.
% 
% IMPORTANT! If you resize the window, you must call TIMEPLT again
% to make the labels look right and have the vector sticks have the
% right orientation.  Matlab will resize the window when you PRINT,
% so you need to ensure that the aspect ratio of the ''paperposition'' 
% is the same as what the aspect ratio of the figure ''position'' is.
% A simple way to do this is to use the FIXPAPER command (included)
% in this distribution) before you print, which will ensure that what
% you get on paper is what you see on the screen!
% 
% -Rich Signell
