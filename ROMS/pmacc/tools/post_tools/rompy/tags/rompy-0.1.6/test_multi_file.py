#!/usr/bin/env python
import glob
import datetime as dt

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import ticker
from matplotlib.colors import Normalize, ListedColormap, LinearSegmentedColormap, hsv_to_rgb
from matplotlib.cm import ScalarMappable

from rompy.extract_from_series import extract_from_series
from rompy import utils
from rompy.plot_utils import banas_hsv_cm, plot_time_series_profile

def my_formatter(x,pos=None):
	return dt.datetime.utcfromtimestamp(x).strftime('%Y-%m-%d %H:%M')



#---------------------------------------------------------------------------------------------------
filelist = glob.glob('ocean_his*.nc')

#x,y = utils.station_to_lat_lon([20])
x = [-122.3904]
y = [47.6108]
varname='salt'
clim_map = {
	'temp':[0,8,20,20],
	'salt':[0,21,33,33],
	'U':[-2,-2,2,2]
	}
#(data, coords) = rompy.extract('ocean_his_1000.nc', varname='salt', extraction_type='profile', x=x, y=y)
(data, time,z) = extract_from_series(filelist,varname=varname,extraction_type='profile',x=x,y=y,freq=dt.timedelta(seconds=1800))

plot_time_series_profile(time,z,data,clim=clim_map[varname],varname=varname)

#print(len(time),len(time[0]))
#print(data)
#print(z)

# 
# fontsize = 8
# 
# clim = [0,21,33,33]
# cmap = banas_hsv_cm(clim[0],clim[1],clim[2],clim[3])
# norm = Normalize(vmin=clim[0],vmax=clim[-1],clip=False)
# sm = ScalarMappable(norm=norm,cmap=cmap)
# sm.set_clim(vmin=clim[0],vmax=clim[-1])
# sm.set_array(np.array([0]))
# 
# 
# fig = Figure(facecolor='white')
# ax1 = fig.add_axes([0.1, 0.55, 0.75, 0.32])
# ax2 = fig.add_axes([0.1, 0.18, 0.75, 0.32])
# cax = fig.add_axes([0.9, 0.1, 0.02, 0.8],frameon=False)
# 
# #ax1 = fig.add_subplot(221)
# #ax1.set_position([0.1, 0.55, 0.75, 0.4])
# #print(ax1)
# #ax2 = fig.add_subplot(223)
# #ax2.set_position([0.1, 0.1, 0.75, 0.4])
# #cax = fig.add_subplot(144,frameon=False)
# #cax.set_position([0.9, 0.1, 0.02, 0.8])
# #fig.subplots_adjust(0,0,1,1)
# 
# my_plot11 = ax1.contourf(time,z,data,norm=norm,cmap=cmap)
# my_plot21 = ax2.contourf(time,z,data,norm=norm,cmap=cmap)
# 
# my_colorbar = fig.colorbar(sm,cax=cax)
# 
# ax1.set_ylim(-20,2)
# ax1.set_xlim(time[0][0],time[-1][-1])
# for yticklabel in ax1.get_yticklabels():
# 		yticklabel.set_fontsize(fontsize)
# ax1.set_xticklabels('')
# 
# ax2.set_xlim(time[0][0],time[-1][-1])
# ax2.set_ylim(np.min(z[0,:]),np.max(z[-1,:]))
# for yticklabel in ax2.get_yticklabels():
# 		yticklabel.set_fontsize(fontsize)
# locs = ax2.get_xticks()
# new_labels = []
# ax2.xaxis.set_major_formatter(ticker.FuncFormatter(my_formatter))
# for label in ax2.get_xticklabels():
# 	label.set_ha('right')
# 	label.set_rotation(30)
# #fig.autofmt_xdate()
# ax1.set_title('Salinity Over Time at a Point')
# FigureCanvas(fig).print_png('/Users/lederer/tmp/rompy.time_series_profile_orig.png')
