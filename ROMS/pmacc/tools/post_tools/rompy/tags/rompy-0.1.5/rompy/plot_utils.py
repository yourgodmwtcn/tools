import datetime as dt
import time

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colors import Normalize, ListedColormap, LinearSegmentedColormap, hsv_to_rgb
from matplotlib.cm import ScalarMappable
from matplotlib import ticker

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

import utils

__version__ = '0.1.5'

def time_series_formatter(x,pos=None):
	return dt.datetime.fromtimestamp(x).strftime('%Y-%m-%d %H:%MZ')

def map_varname(v):
	mapping = {
		'temp':'Temperature',
		'salt':'Salinity',
		'U':'Velocity',
		}
	return mapping[v]

def red_blue_cm():
	cdict = {'red':		[(0.0, 0.0, 0.0),
						(0.5,1.0,1.0),
						(1.0, 0.83, 0.83)],
			
			'green':	[(0.0, 0.34, 0.34),
						(0.5, 1.0, 1.0),
						(1.0, 0.0, 0.0)],
			
			'blue':		[(0.0, 0.75, 0.75),
						(0.5, 1.0, 1.0),
						(1.0, 0.0, 0.0)]
			}
	return LinearSegmentedColormap('red_blue_cm',cdict,N=256)
#	return ListedColormap(['b','w','r'],name='red_blue',N=None)

def banas_cm(a,b,c,d):
	norm = Normalize(vmin=a,vmax=d,clip=False)
	cdict = {'red':[],'green':[],'blue':[]}
	
	if not a==b:
		# add dark blue
		cdict['red'].append((0., 0., 0.))
		cdict['green'].append((0., 0., 0.))
		cdict['blue'].append((0., 0., 0.25))
		
		# add blue
		cdict['red'].append((norm(b), 0., 0.))
		cdict['green'].append((norm(b), 0., 0.))
		cdict['blue'].append((norm(b), 1.0, 1.0))
	else:
		cdict['red'].append((0., 0., 0.))
		cdict['green'].append((0., 0., 0.))
		cdict['blue'].append((0., 0., 1.0))

	# add green between blue and yellow
	cdict['red'].append((norm(b + (c-b)/4.0), 0., 0.))
	cdict['green'].append((norm(b + (c-b)/4.0), 1.0, 1.0))
	cdict['blue'].append((norm(b + (c-b)/4.0), 0., 0.))

	# add yellow in the middle
	cdict['red'].append((norm((b+c)/2.0), 1.0, 1.0))
	cdict['green'].append((norm((b+c)/2.0), 1.0, 1.0))
	cdict['blue'].append((norm((b+c)/2.0), 0., 0.))

	if not c==d:
		# add red
		cdict['red'].append((norm(c), 1.0, 1.0))
		cdict['green'].append((norm(c), 0., 0.))
		cdict['blue'].append((norm(c), 0., 0.))
		
		# add dark red
		cdict['red'].append((1.0, 0.25, 0.25))
		cdict['green'].append((1.0, 0., 0.))
		cdict['blue'].append((1.0, 0., 0.))
	else:
		cdict['red'].append((1.0, 1.0, 1.))
		cdict['green'].append((1.0, 0., 0.))
		cdict['blue'].append((1.0, 0., 0.))
	
	return LinearSegmentedColormap('banas_cm',cdict,N=100)

def banas_hsv_cm(a,b,c,d,N=100):
	norm = Normalize(vmin=a,vmax=d,clip=False)
	cdict = {'red':[],'green':[],'blue':[]}
	if N >= 100:
		n = N
	else:
		n = 100
	aa = norm(a) # 0.0
	bb = norm(b)
	cc = norm(c)
	yy = 0.5*(bb+cc) # yellow is half way between blue and red
	dd = norm(d) # 1.0
	center_value = 0.87
	end_value = 0.65
	tail_end_value = 0.3
	
	blue_hue = 0.55
	yellow_hue = 1./6.
	red_hue = 0.04
	green_hue = 1./3.
	
	gg = ((green_hue - blue_hue)/(yellow_hue - blue_hue))*(yy-bb) + bb
	green_desaturation_width = 0.67
	green_desaturation_amount = 0.5
	
	ii = np.linspace(0.,1.,n)
	hue = np.zeros(ii.shape)
	sat = np.ones(ii.shape)
	val = np.zeros(ii.shape)
	hsv = np.zeros((1,n,3))
	
	val_scaler = -(center_value - end_value)/((cc-yy)*(cc-yy))
	hue_scaler = -(blue_hue - yellow_hue)/((yy-bb)*(yy-bb))
	
	for i in range(len(ii)):
		if ii[i] < bb: # if true then aa is less than bb
			#hue[i] = blue_hue
			hsv[0,i,0] = blue_hue
			#val[i] = tail_end_value*(1 - (ii[i]-aa)/(bb-aa) ) + end_value*( (ii[i]-aa)/(bb-aa) )
			hsv[0,i,2] = tail_end_value*(1 - (ii[i]-aa)/(bb-aa) ) + end_value*( (ii[i]-aa)/(bb-aa) )
		elif ii[i] <= yy:
			#hsv[0,i,0] = blue_hue*(1 - (ii[i]-bb)/(yy-bb) ) + yellow_hue*( (ii[i]-bb)/(yy-bb) )
			hsv[0,i,0] = hue_scaler*(ii[i] -2*bb + yy)*(ii[i] - yy)+yellow_hue
			hsv[0,i,2] = end_value*(1 - (ii[i]-bb)/(yy-bb) ) + center_value*( (ii[i]-bb)/(yy-bb) )
		elif ii[i] <= cc:
			hsv[0,i,0] = yellow_hue*(1 - (ii[i]-yy)/(cc-yy) ) + red_hue*( (ii[i]-yy)/(cc-yy) )
			#hsv[0,i,2] = center_value*(1 - (ii[i]-yy)/(cc-yy) ) + end_value*( (ii[i]-yy)/(cc-yy) )
			hsv[0,i,2] = val_scaler*(ii[i] -2*yy + cc)*(ii[i] - cc)+end_value
		elif ii[i] <= dd:
			hsv[0,i,0] = red_hue
			hsv[0,i,2] = end_value*(1 - (ii[i]-cc)/(dd-cc) ) + tail_end_value*( (ii[i]-cc)/(dd-cc) )
		hsv[0,i,1] = 1.0 - green_desaturation_amount * np.exp(-np.power(3.0*(ii[i]-gg)/((cc-bb)*green_desaturation_width),2.0))
	

#	plt.plot(np.linspace(a,d,n),hsv[0,:,0],'r',np.linspace(a,d,n),hsv[0,:,1],'g',np.linspace(a,d,n),hsv[0,:,2],'b')
#	plt.show()
	
	
	rgb = hsv_to_rgb(hsv)
	cdict['red'].append((0.,0.,rgb[0,0,0]))
	cdict['green'].append((0.,0.,rgb[0,0,1]))
	cdict['blue'].append((0.,0.,rgb[0,0,2]))
	
	for j in range(len(ii)-2):
		i = j+1
		cdict['red'].append((ii[i],rgb[0,i,0],rgb[0,i+1,0]))
		cdict['green'].append((ii[i],rgb[0,i,1],rgb[0,i+1,1]))
		cdict['blue'].append((ii[i],rgb[0,i,2],rgb[0,i+1,2]))

	cdict['red'].append((1.0,rgb[0,-1,0],rgb[0,-1,0]))
	cdict['green'].append((1.0,rgb[0,-1,1],rgb[0,-1,1]))
	cdict['blue'].append((1.0,rgb[0,-1,2],rgb[0,-1,2]))
	
	return LinearSegmentedColormap('banas_cm',cdict,N=N)

def make_cmap_sm_norm(d=None,clim=None,cmap=None):
	if cmap == 'red_blue':
		cmap = red_blue_cm()
	if cmap == 'banas_cm':
		if clim==None:
			cmap = banas_cm(np.min(d[:]),np.min(d[:]),np.max(d[:]),np.max(d[:]))
		elif len(clim) == 2:
			cmap = banas_cm(clim[0],clim[0],clim[1],clim[1])
		elif len(clim) == 4:
			cmap = banas_cm(clim[0],clim[1],clim[2],clim[3])
	elif cmap == 'banas_hsv_cm':
		if clim==None:
			cmap = banas_hsv_cm(np.min(d[:]),np.min(d[:]),np.max(d[:]),np.max(d[:]))
		elif len(clim) == 2:
			cmap = banas_hsv_cm(clim[0],clim[0],clim[1],clim[1])
		elif len(clim) == 4:
			cmap = banas_hsv_cm(clim[0],clim[1],clim[2],clim[3])
	
	norm = Normalize(vmin=clim[0],vmax=clim[-1],clip=False)
	sm = ScalarMappable(norm=norm,cmap=cmap)
	sm.set_clim(vmin=clim[0],vmax=clim[-1])
	sm.set_array(np.array([0]))

	return cmap,sm,norm

def plot_surface(x,y,data,filename='/Users/lederer/tmp/rompy.tmp.png'):
	print('Making plot')
	fig = Figure(facecolor='white',figsize=(12.0,12.0))
	ax = fig.add_subplot(111)
	
	ax.pcolormesh(x,y,data)
#	ax.contour(x,y,data,20)
	ax.axis('tight')
	ax.set_aspect('equal')
	ax.grid()
	FigureCanvas(fig).print_png(filename)

def plot_map(lon,lat,data,filename='/Users/lederer/tmp/rompy.map.png',resolution='h',clim=None,cmap='banas_hsv_cm',title=None, caxis_label=None):
	fig = Figure(facecolor='white',figsize=(12.0,9.0))
#	ax = fig.add_subplot(111)
	longest_side_size = 24.0
	#ax = fig.add_axes((0.,0.,1.,1.),axisbg='grey')
	
	cmap,sm,norm = make_cmap_sm_norm(d=data,clim=clim,cmap=cmap)
	
	ax1 = fig.add_axes((0.1,0.1,0.4,0.8),axisbg='grey')
	ax2 = fig.add_axes((0.5,0.1,0.4,0.8),axisbg='grey')
	cax = fig.add_axes([0.9, 0.1, 0.02, 0.8],frameon=False)
	lllat = np.min(lat)
	urlat = np.max(lat)
	lllon = np.min(lon)
	urlon = np.max(lon)
	
	# puget sound bounding box
	psbb_lllat = 47.0
	psbb_urlat = 48.5
	psbb_lllon = -123.2
	psbb_urlon = -122.1
	
#	print(lllat,urlat,lllon,urlon)
#	m1 = Basemap(projection='merc',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon,resolution=resolution,ax=ax1)
#	m2 = Basemap(projection='merc',llcrnrlat=psbb_lllat,urcrnrlat=psbb_urlat,llcrnrlon=psbb_lllon,urcrnrlon=psbb_urlon,resolution='f',ax=ax2)
	m1 = Basemap(projection='merc',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon,resolution='c',ax=ax1)
	m2 = Basemap(projection='merc',llcrnrlat=psbb_lllat,urcrnrlat=psbb_urlat,llcrnrlon=psbb_lllon,urcrnrlon=psbb_urlon,resolution='c',ax=ax2)
	x1,y1 = m1(*(lon,lat))
	x2,y2 = m2(*(lon,lat))

# Code to make the map fit snuggly with the png
#print(np.max(x), np.min(x), np.max(y),np.min(y))
#	width = np.max(x) - np.min(x)
#	height = np.max(y) - np.min(y)
	
#	if width >= height:
#		fig.set_size_inches(longest_side_size, (height/width)*longest_side_size)
#	else:
#		fig.set_size_inches((width/height)*longest_side_size, longest_side_size)
#	ax.set_position([0.,0.,1.,1.])




# 	bbox = ax.get_position()
# 	print(bbox.xmin, bbox.xmax, bbox.ymin, bbox.ymax)
# 	
# 	fig.set_size_inches((bbox.xmax - bbox.xmin)*longest_side_size, (bbox.ymax - bbox.ymin)*longest_side_size)
# 	ax.set_position([0.,0.,1.,1.])
# 	bbox = ax.get_position()
# 	print(bbox.xmin, bbox.xmax, bbox.ymin, bbox.ymax)
# 	
# 	
# 	if clim==None:
# 		cmap = banas_hsv_cm(np.min(data[:]),np.min(data[:]),np.max(data[:]),np.max(data[:]))
# 		norm = Normalize(vmin=np.min(data[:]),vmax=np.max(data[:]),clip=False)
# 	elif len(clim) == 2:
# 		cmap = banas_hsv_cm(clim[0],clim[0],clim[1],clim[1],N=20)
# 		norm = Normalize(vmin=clim[0],vmax=clim[-1],clip=False)
# 	elif len(clim) == 4:
# 		cmap = banas_hsv_cm(clim[0],clim[1],clim[2],clim[3])
#		norm = Normalize(vmin=clim[0],vmax=clim[-1],clip=False)
		
	pcm1 = m1.pcolormesh(x1,y1,data,cmap=cmap,norm=norm)
#	m1.drawcoastlines(linewidth=0.5)
	regional_coast_lon, regional_coast_lat = m1(*(utils.get_coastline('regional')))
	m1.plot(regional_coast_lon,regional_coast_lat,'k',linewidth=0.5)
	
	pcm2 = m2.pcolormesh(x2,y2,data,cmap=cmap,norm=norm)
#	m2.drawcoastlines(linewidth=0.5)
	detailed_coast_lon, detailed_coast_lat = m2(*(utils.get_coastline('detailed')))
	m2.plot(detailed_coast_lon,detailed_coast_lat,'k',linewidth=0.5)
	
	my_colorbar = fig.colorbar(sm,cax=cax)
	if not caxis_label == None:
		my_colorbar.set_label(caxis_label)
	
	if not title == None:
		ax1.set_title(title)
	
	FigureCanvas(fig).print_png(filename)

def plot_profile(data,depth,filename='/Users/lederer/tmp/rompy.profile.png'):
	fig = Figure()
	ax = fig.add_subplot(111)
	
	ax.plot(data,depth)
	
	ax.grid()
	
	FigureCanvas(fig).print_png(filename)

def plot_mickett(coords,data,varname='',region='',filename='/Users/lederer/tmp/rompy.mickett.png',n=1,x_axis_style='kilometers',x_axis_offset=0,clim=None,cmap=None,labeled_contour_gap=None):
	fig = Figure(facecolor='white')
	fontsize = 8
	
	cmap,sm,norm =  make_cmap_sm_norm(d=data,clim=clim,cmap=cmap)
	
	ax1 = fig.add_axes([0.1, 0.55, 0.75, 0.4])
	ax2 = fig.add_axes([0.1, 0.1, 0.75, 0.4])
	cax = fig.add_axes([0.9, 0.1, 0.02, 0.8],frameon=False)
	
	x_axis_as_km = utils.coords_to_km(coords)
	station_locations = x_axis_as_km[0:-1:n]
# 	
# 	if not clim == None:
# 		norm = Normalize(vmin=clim[0],vmax=clim[-1],clip=False)
# 		sm = ScalarMappable(norm=norm,cmap=cmap)
# 		sm.set_clim(vmin=clim[0],vmax=clim[-1])
# 		sm.set_array(np.array([0]))
# 	else:
# 		norm = None	
# 	
	my_plot11 = ax1.contourf(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,100,norm=norm,cmap=cmap)
	my_plot12 = ax1.contour(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,100,linewidths=1,linestyle=None,norm=norm,cmap=cmap)
	
	if labeled_contour_gap is not None:
		if int(labeled_contour_gap) == labeled_contour_gap:
			contour_label_fmt = '%d'
		else:
			contour_label_fmt = '%1.2f'
		
		solid_contours = np.arange(clim[0],clim[-1],labeled_contour_gap)
#		ax1_xlim = ax1.get_xlim()
		my_plot13 = ax1.contour(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,solid_contours,colors='k',linewidths=0.5)
		ax1.clabel(my_plot13,inline=True,fmt=contour_label_fmt,fontsize=fontsize)
#		ax1.set_xlim(ax1_xlim)

	my_plot14 = ax1.plot(station_locations, 1.5*np.ones(len(station_locations)),'v',color='grey')
	
	ax1.fill_between(x_axis_as_km,coords['zm'][0,:],ax1.get_ylim()[0],color='grey')
#	ax1.set_ylim((-20,ax1.get_ylim()[1]))
	ax1.set_ylim((-20,2))
	ax1.set_xlim((0,x_axis_as_km[-1]))
	for yticklabel in ax1.get_yticklabels():
		yticklabel.set_fontsize(fontsize)
	
	
	
	my_plot21 = ax2.contourf(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,100,norm=norm,cmap=cmap)
	my_plot22 = ax2.contour(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,100,linewidths=1,linestyle=None,norm=norm,cmap=cmap)

	if labeled_contour_gap is not None:
		
#		ax2_xlim = ax2.get_xlim()
		my_plot23 = ax2.contour(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,solid_contours,colors='k',linewidths=0.5)
		ax2.clabel(my_plot23,inline=True,fmt=contour_label_fmt,fontsize=fontsize)
#		ax2.set_xlim = ax2_xlim
#	print(ax2.get_ylim())
#	ax2.fill_between(x_axis_as_km,coords['zm'][0,:],ax2.get_ylim()[0],color='grey')
	ax2.fill_between(x_axis_as_km,coords['zm'][0,:],-1000.0,color='grey')
#	ax2.set_ylim(ax2.get_ylim()[0],2)
#	print(ax2.get_ylim())
	ax2.set_ylim((np.min(coords['zm'][:])-20.0),2)
#	print(ax2.get_ylim())
	
	ax2.set_xlim((0,x_axis_as_km[-1]))
	for yticklabel in ax2.get_yticklabels():
		yticklabel.set_fontsize(fontsize)
# 	if clim == None:
# 		sm = my_plot11
	
	my_colorbar = fig.colorbar(sm,cax=cax)
	if labeled_contour_gap is not None:
		my_colorbar.add_lines(my_plot23)
	
	ax1.set_title('%s %s from a ROMS run' % (region,varname))
	ax1.set_ylabel('depth in meters',position=(0.05,0))
#	ax1.set_xticks(10*np.arange(x_axis_as_km[-1]/10))
	ax1.set_xticks(station_locations)
	ax1.set_xticklabels('')
	
	if x_axis_style == 'kilometers' or x_axis_style == 'kilometer':
			#tick_list = x_axis_as_km[::n]
			#ax2.set_xticks(tick_list)
			#ax2.set_xticklabels([int(tick) for tick in tick_list],size=fontsize)
			td = 10 #tick_distance
			
			ax2.set_xticks(td*np.arange(x_axis_as_km[-1]/td) + (x_axis_offset % td))
			ax2.set_xticklabels([int(num) for num in np.arange(-int(x_axis_offset - x_axis_offset % td),x_axis_as_km[-1],td)])
			for xticklabel in ax2.get_xticklabels():
				xticklabel.set_fontsize(fontsize)
			ax2.set_xlabel('Kilometers')

	elif x_axis_style == 'stations' or x_axis_style == 'station':
		if region == 'Hood Canal':
			tick_list = x_axis_as_km[::n]
			ax2.set_xticks(tick_list)
			ax2.set_xticklabels(utils.hood_canal_station_list(),size=fontsize)
			ax2.set_xlabel('Station ID')
			
		elif region == 'Main Basin':
			tick_list = x_axis_as_km[::n]
			ax2.set_xticks(tick_list)
			ax2.set_xticklabels(utils.main_basin_station_list(),size=fontsize)
	 		ax2.set_xlabel('Station ID') 			
 	else:
 		ax2.set_xticks(x_axis_as_km)
 		ax2.set_xticklabels('')
		ax2.set_xlabel('Kilometers')
	
	FigureCanvas(fig).print_png(filename)

def plot_time_series_profile(t,z,d,filename='/Users/lederer/tmp/rompy.time_series_profile.png',clim=None,cmap='banas_hsv_cm',varname=None, title=None, caxis_label=None):
	
	fontsize = 8
	
	cmap,sm,norm = make_cmap_sm_norm(d=d,clim=clim,cmap=cmap)
	
	
	fig = Figure(facecolor='white')
	ax1 = fig.add_axes([0.1, 0.55, 0.75, 0.32])
	ax2 = fig.add_axes([0.1, 0.18, 0.75, 0.32])
	cax = fig.add_axes([0.9, 0.18, 0.02, 0.69],frameon=False)
	
	my_plot11 = ax1.contourf(t,z,d,100,norm=norm,cmap=cmap)
	my_plot12 = ax1.contour(t,z,d,100,linewidths=1,linestyle=None,norm=norm,cmap=cmap)
	my_plot21 = ax2.contourf(t,z,d,100,norm=norm,cmap=cmap)
	my_plot22 = ax2.contour(t,z,d,100,linewidths=1,linestyle=None,norm=norm,cmap=cmap)

	my_colorbar = fig.colorbar(sm,cax=cax)
	if not caxis_label == None:
		my_colorbar.set_label(caxis_label)
	
	ax1.set_ylim(-20,2)
	ax1.set_xlim(t[0][0],t[-1][-1])
	
	# lets pick some x ticks that aren't stupid
	xmin_dt = dt.datetime.fromtimestamp(t[0][0])
	xmax_dt = dt.datetime.fromtimestamp(t[-1][-1])
	time_window = xmax_dt -xmin_dt
	if (time_window) < dt.timedelta(hours=48):
		date_list = []
		next_time = xmax_dt- dt.timedelta(seconds = xmax_dt.minute*60 + xmax_dt.second)
		while next_time >= xmin_dt:
			date_list.append(next_time)
			next_time = next_time - dt.timedelta(hours=6)
	
	elif (time_window) < dt.timedelta(days=8):
		date_list = []
		next_time = xmax_dt - dt.timedelta(seconds = (xmax_dt.hour*60 + xmax_dt.minute)*60 + xmax_dt.second)
		while next_time >= xmin_dt:
			date_list.append(next_time)
			next_time = next_time - dt.timedelta(days=1)
	
	elif (time_window) < dt.timedelta(days=50):
		date_list = []
		next_time = xmax_dt - dt.timedelta(seconds = (xmax_dt.hour*60 + xmax_dt.minute)*60 + xmax_dt.second)
		while next_time >= xmin_dt:
			date_list.append(next_time)
			next_time = next_time - dt.timedelta(days=7)
	else :
		date_list = [xmin_dt, xmax_dt]
	x_tick_list = []
	for date in date_list:
		x_tick_list.append(time.mktime(date.timetuple()))
	ax2.xaxis.set_major_locator(ticker.FixedLocator(x_tick_list))
	
	for yticklabel in ax1.get_yticklabels():
			yticklabel.set_fontsize(fontsize)
	ax1.set_xticklabels('')
	
	ax2.set_xlim(t[0][0],t[-1][-1])
	ax2.set_ylim(np.min(z[0,:]),np.max(z[-1,:]))
	for yticklabel in ax2.get_yticklabels():
			yticklabel.set_fontsize(fontsize)
	locs = ax2.get_xticks()
	new_labels = []
	ax2.xaxis.set_major_formatter(ticker.FuncFormatter(time_series_formatter))
	for label in ax2.get_xticklabels():
		label.set_ha('right')
		label.set_rotation(30)
	
	if title == None or title == '':
		ax1.set_title('%s Over Time at a Point'% map_varname(varname))
	else:
		ax1.set_title(title)
	
	FigureCanvas(fig).print_png(filename)
	
def plot_parker(coords,data,varname='',title=None,region='',filename='/Users/lederer/tmp/rompy.mickett.png',n=1,x_axis_style='kilometers',resolution='i',x_axis_offset=0,clim=None,cmap=None,labeled_contour_gap=None, caxis_label=None):
	fig = Figure(facecolor='white',figsize=(12.0,9.0))
	fontsize = 8
	
	cmap,sm,norm =  make_cmap_sm_norm(d=data,clim=clim,cmap=cmap)
	
	ax1 = fig.add_axes([0.1, 0.55, 0.65, 0.4]) # top 20 meters
	ax2 = fig.add_axes([0.1, 0.1, 0.65, 0.4]) # full column
	ax3 = fig.add_axes([0.7, 0.55, 0.3, 0.4],axis_bgcolor='white')#'#298FAF') # map of domain containing curtain
	cax = fig.add_axes([0.84, 0.1, 0.02, 0.4],frameon=False) # subplot for the color axis
	
	x_axis_as_km = utils.coords_to_km(coords)
	station_locations = x_axis_as_km[0:-1:n]

	my_plot11 = ax1.contourf(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,100,norm=norm,cmap=cmap)
	my_plot12 = ax1.contour(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,100,linewidths=1,linestyle=None,norm=norm,cmap=cmap)
	
	if labeled_contour_gap is not None:
		if int(labeled_contour_gap) == labeled_contour_gap:
			contour_label_fmt = '%d'
		else:
			contour_label_fmt = '%1.2f'
		
		solid_contours = np.arange(clim[0],clim[-1],labeled_contour_gap)

		my_plot13 = ax1.contour(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,solid_contours,colors='k',linewidths=0.5)
		ax1.clabel(my_plot13,inline=True,fmt=contour_label_fmt,fontsize=fontsize)


	my_plot14 = ax1.plot(station_locations, 1.5*np.ones(len(station_locations)),'v',color='grey')
	
	ax1.fill_between(x_axis_as_km,coords['zm'][0,:],ax1.get_ylim()[0],color='grey')

	ax1.set_ylim((-20,2))
	ax1.set_xlim((0,x_axis_as_km[-1]))
	for yticklabel in ax1.get_yticklabels():
		yticklabel.set_fontsize(fontsize)
	
	
	
	my_plot21 = ax2.contourf(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,100,norm=norm,cmap=cmap)
	my_plot22 = ax2.contour(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,100,linewidths=1,linestyle=None,norm=norm,cmap=cmap)

	if labeled_contour_gap is not None:
		my_plot23 = ax2.contour(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,solid_contours,colors='k',linewidths=0.5)
		ax2.clabel(my_plot23,inline=True,fmt=contour_label_fmt,fontsize=fontsize)

	ax2.fill_between(x_axis_as_km,coords['zm'][0,:],-1000.0,color='grey')
	ax2.set_ylim((np.min(coords['zm'][:])-20.0),2)
	
	ax2.set_xlim((0,x_axis_as_km[-1]))
	for yticklabel in ax2.get_yticklabels():
		yticklabel.set_fontsize(fontsize)
	
	my_colorbar = fig.colorbar(sm,cax=cax)
	if not caxis_label == None:
		my_colorbar.set_label(caxis_label)
	
	if labeled_contour_gap is not None:
		my_colorbar.add_lines(my_plot23)
	
	if title==None:
		ax1.set_title('%s %s from a ROMS run' % (region,varname))
	else:
		ax1.set_title(title)
	
	ax1.set_ylabel('depth in meters',position=(0.05,0))
	ax1.set_xticks(station_locations)
	ax1.set_xticklabels('')
	
	if x_axis_style == 'kilometers' or x_axis_style == 'kilometer':
			td = 10 #tick_distance
			
			left_most_tick_label = -x_axis_offset + (x_axis_offset % td)
			left_most_tick = left_most_tick_label + x_axis_offset
						
			ax2.set_xticks(np.arange(left_most_tick,x_axis_as_km[-1],td))
			ax2.set_xticklabels([int(num) for num in np.arange(left_most_tick_label, x_axis_as_km[-1],td)])
			
#			ax2.set_xticks(td*np.arange(x_axis_as_km[-1]/td) + (x_axis_offset % td))
#			ax2.set_xticklabels([int(num) for num in np.arange(-int(x_axis_offset - x_axis_offset % td),x_axis_as_km[-1],td)])
			for xticklabel in ax2.get_xticklabels():
				xticklabel.set_fontsize(fontsize)
			ax2.set_xlabel('Kilometers')

	elif x_axis_style == 'stations' or x_axis_style == 'station':
		if region == 'Hood Canal':
			tick_list = x_axis_as_km[::n]
			ax2.set_xticks(tick_list)
			ax2.set_xticklabels(utils.hood_canal_station_list(),size=fontsize)
			ax2.set_xlabel('Station ID')
			
		elif region == 'Main Basin':
			tick_list = x_axis_as_km[::n]
			ax2.set_xticks(tick_list)
			ax2.set_xticklabels(utils.main_basin_station_list(),size=fontsize)
	 		ax2.set_xlabel('Station ID') 			
 	else:
 		ax2.set_xticks(x_axis_as_km)
 		ax2.set_xticklabels('')
		ax2.set_xlabel('Kilometers')
	
	# make map in the top right corner
	# these lat lon values are derived from the curtain defined for the plot
#	lllat = np.min(coords['ym'])
#	urlat = np.max(coords['ym'])
#	lllon = np.min(coords['xm'])
#	urlon = np.max(coords['xm'])
	
#	lat lon values for the inset map show a close-up of the Puget Sound
	lllat = 47.0
	urlat = 48.5
	lllon = -123.3
	urlon = -122.2
	
	
#	m = Basemap(projection='merc',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon,resolution=resolution,ax=ax3)
	m = Basemap(projection='merc',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon,resolution='c',ax=ax3)
	x,y = m(*(coords['xm'],coords['ym']))
#	pcm = m.plot(x,y,'r')
	coast_lon, coast_lat = m(*(utils.get_coastline('detailed')))
#	print(coast_lat)
#	print(coast_lon)
#	m.drawcoastlines(linewidth=0.5)
#	m.fillcontinents(color='#ECECEC')
	m.plot(coast_lon,coast_lat,'k',linewidth=0.5)
	pcm1 = m.plot(x,y,'r',linewidth=0.5)
	pcm2 = m.plot(x[0:-1:n],y[0:-1:n],'.k')
	
	FigureCanvas(fig).print_png(filename)
	