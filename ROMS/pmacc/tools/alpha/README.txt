*************** README for the "tools" **********************
Parker MacCready, University of Washington, 2012-

GENERAL INSTRUCTIONS:

See: faculty.washington.edu/pmacc/tools/modeling_tools.html

NOTES (most recent on top):

1/16/2013 I recently recoded my mooring extraction code, changing the overall organization.  Now all mooring extraction and plotting are done from the folder /tools/moor (extracted files go to /tools_output/moor_out).  This is all assuming you use my tools setup. You would need also alpha, pandora, shared, and port_tools. The extraction is run by get_moor.m (calls Z_get_moor.m), and plotting is done by moor_plot.m, which can call a number of different plotting routines. Also I re-coded Z_get_moor.m so that it can accept a vector of lat, lon locations (provided by get_moor.m) instead of just a single point.  Then in the code it extracts and stores the info for all mooring locations while opening each history file just once.  The idea was to try to speed up large extractions for minimizing the time for opening and closing all the netcdf history files.  It stores the result in a structure M (just like mod_moor in earlier versions but easier to type!), and M is a structure array.  So for example to get the zeta time series for the 3rd mooring in the list you would use M(3).zeta. I tested it and it seems to work, although the speed was not as fabulous as I  had hoped.  On euclid (a NERSC machine in Berkeley) it did three mooring locations for a full year-of-hours PNWTOX 40 layer run (no dye or bio) in 3.7 hours, and it did one mooring in just over an hour, so the scaling seems OK.  On skua working on salish_2006_4 it took 11 hours!  I'm not sure why skua is so slow.

1/14/2013 Created the separate folder tools/moor/ do deal with all mooring extraction and plotting code.  This includes putting the main extraction function tools/pandora/Z_functions/Z_get_moor.m into tools/moor/Z_functions/Z_get_moor.m. Results go into tools_output/moor_out/.

1/14/2013 Replaced the folder t_tide_v1.1/ with the newer t_tide_v1/ downloaded today from Rich Pawlowicz's page: http://www.eos.ubc.ca/~rich/.  This gets rid of a number of error messages resulting from calls to out-of-date functions.

1/7/2013 I recently changed the rtools architecture to include a new directory tools/rtools_user/.  This is meant to be the place where a user of rtools would keep their own personal grid and run definitions, so that they will not be overwritten when downloading a newer version of rtools/.

12/19/2012 I created Z_s2z.m (based on Sarah Giddings Z_s2z_mat2.m), and it sits in pandora/Z_functions. I then went through all the code in tools/ and changed the function call.  I also edited Z_get_basic_info so that it gets the variable S.Vtransform that Z_s2z.m requires.  Then I got rid of all other instances of this code: pandora/Z_functions/Z_s2z_mat.m, rtools/Z_utils/Z_s2z_mat.m and rtools/Z_utils/Z_s2z_mat2.m.  There is still the code post_tools/utility/roms_z.m which does the same thing, 3-4 times faster, but which will give correct results only for hc=0 and Vtransform=1.  It is still used in many places.  The results of the two were identical when tested on a salish_2006_4 history file.
