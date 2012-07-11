% -------------------------------------------------------------------------
%                       RANDOM
% ------------------------------------------------------------------------- 
% avg1                      - calculates average of all sucessive 2 point intervals
% beautify_session          - modifies various figure/axes/line/text object
%                             properties for the session.
% clr                       - clears all variables and closes all figures
%
% sdisp                     - Pretty print simplified symbolic expression
%
% calc_error                - returns maximum percentage error (in decimal) 
% 
% check_gap                 - 
% fill_data                 - replaces NaN's with linearly interpolated values
% fill_gap                  -
% find_gap                  -
%
% find_file                 - finds valid filenames using regexp();
%
% simple_ls                 - simple least squares by right division
% taper_ls                  - smoothed & tapered least squares
% svd_fit                   - svd solution of system
%
% addnan                    - replaces values > input with NaN.
% fillnan                   - replaces values == input with NaN. More reasonable 
% repnan                    - Replaces NaN's with specified value
%
% gen_ser                   - generates a series 1:size(variable,dimension)
% nan_detrend               - Removes columnwise mean from each column of var.
% revz                      - reverses ydir
% stat                      - returns variable statistics. Accepts string input.
% txp                       - calculates 10^(). Accepts string input
% vecrotate                 - Rotates vectors by angle theta in degrees.
%
% -------------------------------------------------------------------------
%                       OCEANOGRAPHY
% -------------------------------------------------------------------------
% bfrq                      - calculates buoyancy frequency
% vertmode                  - calculates vertical modes
% pv                        - calculates pv
%
% roms_movie                - displays movie of variable in ROMS output file
% roms_diffmovie            - displays movie of variable - variable at
%                             first "time" instant.
% roms_pv                   - calculates Ertel PV from ROMS output file
% roms_energy               - 
%
% -------------------------------------------------------------------------
%                       TIME SERIES STUFF
% -------------------------------------------------------------------------
% aliasfreq                 - returns apparent frequency given a known frequency
% confchi2                  - Confidence interval for chi-squared variate
% conft                     - Confidence interval for t-distributed variate
% dcdetide                  - removes the main 4 tidal constituents
% dcconv                    - calculates discrete convolution (1D only). replicates MATLAB's conv. use that instead.
% dccoher                   - Calculates coherence amplitude and phase.
% dcfft                     - calculates fft and returns power, frequency, FFT coeff's
% dchist2                   - plots 2d histogram (BETTER THAN BELOW)
% dchist2d                  - calls hist2d (in misc) with appropriate arguments and plots result.
% dcpgram                   - plots periodogram of input. does fft inside and marks top 'x' peaks
% dcpsd                     - computes and plots power spectrum density
% dcstructfn                - computes and plots structure function
% chkparsvl                 - takes fft coeffs and data and checks whether parseval's  theorem is satisfied. 
%                           - dcfft calls this by default.
%
% -------------------------------------------------------------------------
%                           PLOT / COLOR STUFF
% -------------------------------------------------------------------------
% animate                   - Shows movie of input data
% animate_quiver            - Shows movie of vectors
% beautify                  - call after making plot to make kickass
% datex                     - calls datetick on x-axis
% disp_plot                 - plots time series with depth displaced by value of depth
% liney                     - plots and labels horizontal line at given y
% loglinex(x,factor)        - Plots and labels vertical line at a given x and a factor 
%                             to scale the labelled 1./x value by.%
% fix_subplot2x2            - Tries to reduce extra whitespace in 2x2 subplots. Params
%                             might need to be set manually
%
%
% cbrewer                   - color brewer
% haxby                     - bathymetric colormap
% lbmap                     - Light & Bartlein Eos article
% cptcmap                   - GMT colormaps for matlab
% diverging_map             - creates divergent colormap based on http://www.sandia.gov/~kmorel/documents/ColorMaps/
% more colormaps            - http://soliton.vm.bytemark.co.uk/pub/cpt-city/index.html
% distinguishable_colors    - cycles best possible colors
% laprint                   - latex print stuff + psfrag stuff
% cm_and_cb_utilities       - colormap and colorbar  = cbfreeze, cblabel etc.
% freezeColors              - freeze colormap etc.
%
%
% quiverc                   - color coded quiver plot
% quiverscale               - plots a reference arrow for quiver plots
% fixepsbbox                - fixes .eps bounding box to reduce extra whitespace
% suplabel                  - super label for subplots
% plotyyy                   - 3 y axes
% legendflex                - flexible legends
% scrollsubplot             - subplots with scrollbars.
% magnifyOnFigure           - is what is says
% pptfigure                 - 
% export_fig                - good figure saving
%
% -------------------------------------------------------------------------
%                           MISC
% ------------------------------------------------------------------------- 
% cprintf                   - colored text output
%
% find_approx               - find approximate equality
% hist2d                    - plot 2d histogram
% matlab2latex              - is what it says
% run_avg2                  - Tom's running average code
% llpa                      - Ken's LLPA low pass filter for detiding.
% vgrid                     - Ken's ROMS vertical grid code.
% spice                     - calculate spiciness
% intersections             - find intersections of two curves
% timeit                    - function run time
%
% fdep                      - lists dependencies for input .m file
% exportToZip               - Creates a ZIP file containing all dependencies of a function
% v2struct                  - variable to structure and vice versa
%
% deg2utm                   - lat/lon to UTM
% 
% -------------------------------------------------------------------------
%                       TOOLBOXES
% -------------------------------------------------------------------------
% timeutil                  - time / date utilities 
%                             http://home.online.no/~pjacklam/matlab/software/util/index.html
% timeplt                   - 
% movieman                  - Ryan's movie stuff
% ROMS                      - roms_wilkin, arango, ROMS_tools
% tsplot                    - 

help dctools