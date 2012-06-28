#!/usr/bin/env python
import datetime as dt

from rompy import extract_from_series
from rompy.plot_utils import plot_time_series_profile
from rompy.utils import station_to_lat_lon

#x = [-122.3904]
#y = [47.6108]

#x,y = station_to_lat_lon(['admiralty inlet'])
x,y = station_to_lat_lon([4])


dt0 = dt.datetime(2006,06,11,0,0,0)
dt1 = dt.datetime(2006,06,13,0,0,0)

interval = 3600 # seconds

varname='salt'

clim_map = {
	'temp':[0,8,20,20],
	'salt':[0,21,33,33],
	'U':[-2,-2,2,2]
	}

d,t,z = extract_from_series.extract_from_two_datetimes(x,y,dt0,dt1,varname=varname,interval=interval)

#print('z')
#print(z)
#print('time')
#print(t)
#print('data')
#print(d)
plot_time_series_profile(t,z,d,clim=clim_map[varname],varname=varname,filename='/Users/lederer/tmp/rompy.arbitrary_time_series_profile.png')
