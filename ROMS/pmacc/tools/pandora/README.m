% READme.m 8/6/2011 Parker MacCready
%
% about the files and folders in "pandora_code"
%
% This code is meant for postprocessing and plotting of ROMS history files,
% and was developed mainly for the RISE, MoSSea, and PNWTOX runs.  However
% I hope that many of the tools are more general, or are easily adapted to
% other regions.
%
% Quirks:
%
% 1. We generaly store our ROMS output with one time step per history
% file.  This is because our domains are usually huge and so a single time
% snapshot can be very big.  The downside is that forming time series, or
% filtering in time, can be slow.  Our current default is one save per
% hour.
%
% 2. You need various toolboxes:
%    tools/shared/Z_functions (used to simplify plotting code)
%    mexnc (to enable reading NetCDF)
%    snctools (the toolbox of NetCDF commands, like "nc_varget"
%    psvs (Neil Banas' toolbox of ROMS routines, very useful)
%    seawater (equation of state for seawater, e.g. sw_dens.m)
%    t_tide (Rich Pawlowicz' tidal analysis toolbox)
%
% 3. Generally the code should work on raw history files (in OUT/) or on
% tidally-averaged files (in OUT_lp).
%
% ** pan_plot.m **
% This is a driver for the plotting code.  It allows you to interactively
% choose the history file (or sequence of files so you can make a movie),
% and then interactively choose the code to plot it with.  Depending on
% which plotting code you use you may also be prompted to select a mooring
% file.  By minor editing it can also be run with the choices hardwired,
% which may be better when runnning on a remote linux machine.  After it
% makes a plot the user can choose to save it or not (by default as a tiff
% file in ../pandora_data/Salish_figures).
% If you chose a sequence of files then it
% automatically creates a directory inside ../pandora_data/Salish_movies
% and saves a sequence of numbered jpeg images there).
%
% ** mooring_extractor.m **
% This a a preprocessing tool.  It extracts a mooring-like record (e.g.
% u(z,t), salt(z,t), etc.) at a user-specified lat-lon location.  It
% extracts only the values in a little square around the mooring before
% doing the x,y interpolation, and this makes it relatively fast.  But it
% still takes about a second per save for a MoSSea run (salish_2006_3, from
% Sutherlad et al. 2011, JPO), and so it takes over 3 hours for a full
% year.  It has an internal switch allowing it to be run in
% interactive/hardwired mode - and hardwired is the way to do it on a
% remote linix machine).  You have to have at least one extracted mooring
% file in order for the plotting code "basic_mooring.m" to run.
%
% ** zslice_preprocess.m **
% This needs to be run before you can use the plotting code
% "basic_zslice.m"  It creates interpolation matrices to speed up the
% creation of sections at a given depth.  Note that it assumes the surface
% height to be zero when figuring out this interpolation matrix.  Thus the
% same matrix can be used at any time.  Caveat: the resulting field won't
% be exactly at the specified depth, because of tidal oscillations.  In my
% opinion this reference frame is often better than a strict Eulerian
% average, especially in shallower water.
%
% *** PLOTTING CODE (in Salish_plot_code, but run using pan_plot.m) ***
%
% 0. basic.m
% Makes maps of surface fields of salinity and temperature.  Good for a
% quick look.  Has a nice overlay of velocity vectors and a velocity scale
% arrow.
%
% 1. PS_JdF_sect.m
% Makes a map of surface salinity, and then a section of salinity and
% velocity down the thalweg of JdF and PS (uses Z_make_sect.m).  Best to
% use with tidally-averaged (lp) field.
%
% 2. basic_zslice.m
% Makes maps of a property at different depth (uses the results of
% zslice_preprocess.m).
%
% 3. basic_mooring.m
% Makes some maps and then panels of time series (uses the results of
% mooring_extractor.m).
%
% 4. basicPS_Fb_a.m
% Makes maps of depth-integrated buoyancy flux per unit horizontal area,
% over the whole domain and in Puget Sound (hence the prefix "basicPS")
% 
% 5. ...anything alse you make will be available...