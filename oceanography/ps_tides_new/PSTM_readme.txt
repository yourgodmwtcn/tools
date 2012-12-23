This code is a MATLAB version of the fortran code documented in:

Lavelle, J. W., H. O. Mofjeld, E. Lempriere-Doggett, G. A. Cannon,
D. J. Pashinski, E. D. Cokelet, L. Lytle, and S. Gill, 1988:
A multiply-connected channel model of tides and tidal currents in
Puget Sound, Washington and a comparison with updated observations.
NOAA Tech. Memo. ERL PMEL-84,
Pacific Marine Environmental Laboratory, NOAA.

It gives numerical predictions of the section-averaged tidal currents at
many different sections across Puget Sound, as well as the tidal height.

It appears to be accurate, including nodal corrections and being tuned
to observed patterns of surface height, but it comes with NO GUARANTEES
of any sort.

The original MATLAB gui was coded by Dave Winkel (UW-APL) in 2001.

Edited by Parker MacCready (UW-Oceanography) in 2009 to improve
spacing of the gui boxes and include the MATLAB datenum-format
time vector in the saved output.

NOTES:

To run it, just execute PS_tides in MATLAB.  Resize the resulting window as
needed and start making your choices.

It will plot graphical output to the screen, and you can choose to save the output
to a mat-file.

Output velocity is averaged across the chosen section, FLOOD is POSITIVE,
and the units are m s-1.

Output height has units of meters.

Output time zone is always in PACIFIC STANDARD TIME, so be sure to correct
as needed.  It is given in both Julian Days from the start of the chosen year, and
in datenum format.  Note that Julian Days start at 1, not zero.

Cheers,

Parker MacCready 2/26/2009
parker@ocean.washington.edu
(206) 685-9588

Please send me any corrections!