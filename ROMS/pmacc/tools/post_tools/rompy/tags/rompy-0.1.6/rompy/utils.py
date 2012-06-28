import os

import numpy as np
from matplotlib.mlab import griddata
import netCDF4 as nc

__version__ = '0.1.6'

def interp_2d_latlon(lat,lon,data,lati,loni):
	return griddata(lat.reshape(lat.size),lon.reshape(lon.size),data.reshape(data.size),lati,loni)

def interp_2d_xy(x,y,data,xi,yi):
	try:
		di = griddata(x.reshape(x.size),y.reshape(y.size),data.reshape(data.size),xi,yi)
	except TypeError:
		di = np.zeros(xi.size)
		if x.ndim ==2 and y.ndim ==2:
			x_vec = x[0,:]
			y_vec = y[:,0]
		elif x.ndim == 3 and y.ndim == 3:
			x_vec = x[0,0,:]
			y_vec = y[0,:,0]
		else:
			x_vec = x
			y_vec = y
		xl = np.nonzero(x_vec <= xi)[0][-1]
		xh = np.nonzero(xi <= x_vec)[0][0]
		yl = np.nonzero(y_vec <= yi)[0][-1]
		yh = np.nonzero(yi <= y_vec)[0][0]
				
		if not x_vec[xl] == x_vec[xh]:
			xd = (xi-x_vec[xh])/(x_vec[xl]-x_vec[xh])
		else:
			xd = 1.
		if not y_vec[yl] == y_vec[yh]:
			yd = (yi-y_vec[yh])/(y_vec[yl]-y_vec[yh])
		else:
			yd = 1.

		w0 = data[yl,xl]*(1-yd) + data[yh,xl]*yd
		w1 = data[yl,xh]*(1-yd) + data[yh,xh]*yd
	
		di = w0*(1-xd) + w1*xd
	return di

def interp_2d_point(x,y,d,xi,yi):
	if not x.ndim == 1 or not y.ndim == 1:
		raise(TypeError,'interp_2d_from_point needs the x and y to be vectors')
	if not xi.size == 1 or not yi.size == 1:
		print(xi,yi,zi)
		raise(TypeError, 'interp_2d_from_point needs xi and yi to be a single value')
	try:
		xl = np.nonzero(x <= xi)[0][-1]
		xh = np.nonzero(xi <= x)[0][0]
		yl = np.nonzero(y <= yi)[0][-1]
		yh = np.nonzero(yi <= y)[0][0]
	except IndexError, e:
		print('x, xi')
		print(x,xi)
		print('y, yi')
		print(y,yi)
	
	if not x[xl] == x[xh]:
		xd = (xi-x[xh])/(x[xl]-x[xh])
	else:
		xd = 1.
	if not y[yl] == y[yh]:
		yd = (yi-y[yh])/(y[yl]-y[yh])
	else:
		yd = 1.
	
	w0 = d[yl,xl]*(1-yd) + d[yh,xl]*yd
	w1 = d[yl,xh]*(1-yd) + d[yh,xh]*yd
	
	return w0*(1-xd) + w1*xd

def interp_2d_from_list_of_points(x,y,z,d,p_list):
	di = np.zeros(len(p_list),1)
	for i in range(len(p_list)):
		di[i] = interp_2d_point(x,y,d,p_list[i][0],p_list[i][1])
	return di

def interp_2d(x,y,data,xi,yi):
	if x.shape == y.shape and x.shape == data.shape:
		# assume x and y are the same everywhere in their respective dimension
		if x.ndim == 2:
			x_vec = x[0,:]
		else:
			x_vec = x
		if y.ndim == 2:
			y_vec = y[:,0]
		else:
			y_vec = y
		
		# assume xi and yi are vectors
		if xi.ndim == 1 and yi.ndim == 1:
			di = np.zeros((len(yi),len(xi)))
			for i in range(len(xi)):
				for j in range(len(yi)):
					di[j,i] = interp_2d_point(x_vec,y_vec,data,xi[i],yi[j])

		# if xi and yi are not vectors, then lets just do everything on a point by point basis.
		elif xi.shape == yi.shape:
			di = np.zeros(xi.shape)
			for i in range(xi.size):
				di.flat[i] = interp_2d_point(x_vec,y_vec,data,xi.flat[i],yi.flat[i])
		return di
		
	else:#elif (len(x),len(y)) == data.shape:
		print('Do this the other way')

def interp_3d_point(x,y,z,d,xi,yi,zi):
	if not x.ndim == 1 or not y.ndim == 1 or not z.ndim == 1:
		raise(TypeError,'interp_3d_from_point needs the x, y, and z to be vectors')
	if not xi.size == 1 or not yi.size == 1 or not zi.size ==1:
		print(xi,yi,zi)
		raise(TypeError, 'interp_3d_from_point needs xi, yi, and zi to be a single value')
	try:
		xl = np.nonzero(x <= xi)[0][-1]
		xh = np.nonzero(xi <= x)[0][0]
		yl = np.nonzero(y <= yi)[0][-1]
		yh = np.nonzero(yi <= y)[0][0]
		zl = np.nonzero(z <= zi)[0][-1]
		zh = np.nonzero(zi <= z)[0][0]
	except IndexError, e:
		print('x, xi')
		print(x,xi)
		print('y, yi')
		print(y,yi)
		print('z, zi')
		print(z,zi)
#	print((xl,xi, xh),(yl, yi, yh),( zl, zi, zh))
	
	if not x[xl] == x[xh]:
		xd = (xi-x[xh])/(x[xl]-x[xh])
	else:
		xd = 1.
	if not y[yl] == y[yh]:
		yd = (yi-y[yh])/(y[yl]-y[yh])
	else:
		yd = 1.
	if not z[zl] == z[zh]:
		zd = (zi-z[zh])/(z[zl]-z[zh])
	else:
		zd = 1.
	
	i0 = d[zl,yl,xl]*(1-zd) + d[zh,yl,xl]*zd
	i1 = d[zl,yh,xl]*(1-zd) + d[zh,yh,xl]*zd
	
	j0 = d[zl,yl,xh]*(1-zd) + d[zh,yl,xh]*zd
	j1 = d[zl,yh,xh]*(1-zd) + d[zh,yh,xh]*zd
	 
	w0 = i0*(1-yd) + i1*yd
	w1 = j0*(1-yd) + j1*yd
	
	# cludge alert
	if np.abs(w0) > 1.0e35 or np.abs(w1) > 1.0e35:
		return np.nan
	else:
		return w0*(1-xd) + w1*xd

def interp_3d_from_list_of_points(x,y,z,d,p_list):
	di = np.zeros(len(p_list),1)
	for i in range(len(p_list)):
		di[i] = interp_3d_point(x,y,z,d,p_list[i][0],p_list[i][1],p_list[i][2])
	return di

	
def interp_3d(x,y,z,data,xi,yi,zi):
	# we make a lot of assumptions about the incoming data. this is not a universal interpn
	if x.shape == y.shape and x.shape == z.shape and x.shape == data.shape:
		#print('Do this the hard way')
		# assume x, y, and z are the same everywhere in their respective dimension
		if x.ndim == 3:
			x_vec = x[0,0,:]
		else:
			x_vec = x
		if y.ndim == 3:
			y_vec = y[0,:,0]
		else:
			y_vec = y
		if z.ndim == 3:
			z_vec = z[:,0,0]
		else:
			z_vec = z
		
		# assume xi, yi, and zi are vectors
		if xi.ndim == 1 and yi.ndim == 1 and zi.ndim == 1:
			di = np.zeros((len(zi), len(yi),len(xi)))
			for i in range(len(xi)):
				for j in range(len(yi)):
					for k in range(len(zi)):
						di[k,j,i] = interp_3d_point(x_vec,y_vec,z_vec,data,xi[i],yi[j],zi[k])

		# if xi, yi, and zi are not vectors, then lets just do everything on a point by point basis.
		elif xi.shape == yi.shape and xi.shape == zi.shape:
			#print("I'm in the right place")
			di = np.zeros(xi.shape)
			for i in range(xi.size):
				di.flat[i] = interp_3d_point(x_vec,y_vec,z_vec,data,xi.flat[i],yi.flat[i],zi.flat[i])
		return di
		
	elif (len(x),len(y),len(z)) == data.shape:
		print('Do this the other way')


def meshgrid(*xi,**kwargs):
    """
    Return coordinate matrices from one or more coordinate vectors.

    Make N-D coordinate arrays for vectorized evaluations of
    N-D scalar/vector fields over N-D grids, given
    one-dimensional coordinate arrays x1, x2,..., xn.

    Parameters
    ----------
    x1, x2,..., xn : array_like
        1-D arrays representing the coordinates of a grid.
    indexing : 'xy' or 'ij' (optional)
        cartesian ('xy', default) or matrix ('ij') indexing of output
    sparse : True or False (default) (optional)
         If True a sparse grid is returned in order to conserve memory.
    copy : True (default) or False (optional)
        If False a view into the original arrays are returned in order to
        conserve memory

    Returns
    -------
    X1, X2,..., XN : ndarray
        For vectors `x1`, `x2`,..., 'xn' with lengths ``Ni=len(xi)`` ,
        return ``(N1, N2, N3,...Nn)`` shaped arrays if indexing='ij'
        or ``(N2, N1, N3,...Nn)`` shaped arrays if indexing='xy'
        with the elements of `xi` repeated to fill the matrix along
        the first dimension for `x1`, the second for `x2` and so on.

    See Also
    --------
    index_tricks.mgrid : Construct a multi-dimensional "meshgrid"
                     using indexing notation.
    index_tricks.ogrid : Construct an open multi-dimensional "meshgrid"
                     using indexing notation.

    Examples
    --------
    >>> x = np.linspace(0,1,3)   # coordinates along x axis
    >>> y = np.linspace(0,1,2)   # coordinates along y axis
    >>> xv, yv = meshgrid(x,y)   # extend x and y for a 2D xy grid
    >>> xv
    array([[ 0. ,  0.5,  1. ],
           [ 0. ,  0.5,  1. ]])
    >>> yv
    array([[ 0.,  0.,  0.],
           [ 1.,  1.,  1.]])
    >>> xv, yv = meshgrid(x,y, sparse=True)  # make sparse output arrays
    >>> xv
    array([[ 0. ,  0.5,  1. ]])
    >>> yv
    array([[ 0.],
           [ 1.]])

    >>> meshgrid(x,y,sparse=True,indexing='ij')  # change to matrix indexing
    [array([[ 0. ],
           [ 0.5],
           [ 1. ]]), array([[ 0.,  1.]])]
    >>> meshgrid(x,y,indexing='ij')
    [array([[ 0. ,  0. ],
           [ 0.5,  0.5],
           [ 1. ,  1. ]]),
     array([[ 0.,  1.],
           [ 0.,  1.],
           [ 0.,  1.]])]

    >>> meshgrid(0,1,5)  # just a 3D point
    [array([[[0]]]), array([[[1]]]), array([[[5]]])]
    >>> map(np.squeeze,meshgrid(0,1,5))  # just a 3D point
    [array(0), array(1), array(5)]
    >>> meshgrid(3)
    array([3])
    >>> meshgrid(y)      # 1D grid; y is just returned
    array([ 0.,  1.])

    `meshgrid` is very useful to evaluate functions on a grid.

    >>> x = np.arange(-5, 5, 0.1)
    >>> y = np.arange(-5, 5, 0.1)
    >>> xx, yy = meshgrid(x, y, sparse=True)
    >>> z = np.sin(xx**2+yy**2)/(xx**2+yy**2)
    """
    copy = kwargs.get('copy',True)
    args = np.atleast_1d(*xi)
    if not isinstance(args, list):
        if args.size>0:
            return args.copy() if copy else args
        else:
            raise TypeError('meshgrid() take 1 or more arguments (0 given)')

    sparse = kwargs.get('sparse',False)
    indexing = kwargs.get('indexing','xy') # 'ij'


    ndim = len(args)
    s0 = (1,)*ndim
    output = [x.reshape(s0[:i]+(-1,)+s0[i+1::]) for i, x in enumerate(args)]

    shape = [x.size for x in output]

    if indexing == 'xy':
        # switch first and second axis
        output[0].shape = (1,-1) + (1,)*(ndim-2)
        output[1].shape = (-1, 1) + (1,)*(ndim-2)
        shape[0],shape[1] = shape[1],shape[0]

    if sparse:
        if copy:
            return [x.copy() for x in output]
        else:
            return output
    else:
        # Return the full N-D matrix (not only the 1-D vector)
        if copy:
            mult_fact = np.ones(shape,dtype=int)
            return [x*mult_fact for x in output]
        else:
            return np.broadcast_arrays(*output)


def ndgrid(*args,**kwargs):
    """
    Same as calling meshgrid with indexing='ij' (see meshgrid for
    documentation).
    """
    kwargs['indexing'] = 'ij'
    return meshgrid(*args,**kwargs)

def station_to_lat_lon(station_list):
	lon_dict = {7.5: -122.6418, 8.5: -122.5848, 4: -122.542, 6: -122.4932, 7: -122.62009999999999, 8: -122.6053, 9: -122.6673, 10: -122.71980000000001, 11: -123.1375, 12: -123.1142, 13: -123.008, 14: -122.93989999999999, 15: -122.8601, 16: -122.7578, 17: -122.76139999999999, 18: -122.6169, 19: -122.6318, 20: -122.6848, 21: -122.85039999999999, 22: -123.0189, 36.299999999999997: -122.91270574957041, 24: -123.1249, 11.199999999999999: -123.0305, 27: -122.41030000000001, 28: -122.45440000000001, 29: -122.44329999999999, 30: -122.4084, 31: -122.34869999999999, 32: -122.4427, 33: -122.5008, 34: -122.5372204355572, 35: -122.6519213323306, 36: -122.82570322014899, 35.600000000000001: -122.7313463437406, 11.1: -123.06376861727099, 403.10000000000002: -122.8790276962233, 35.5: -122.68581305397331, 'admiralty inlet': -122.7557, 33.5: -122.5621, 34.5: -122.6043449761104, 402.10000000000002: -123.00131004223751, 36.200000000000003: -122.91573023207449, 401: -123.05670000000001, 36.100000000000001: -122.8573269666242, 402: -123.0167366137608, 403: -122.91127210926319}
	lat_dict = {7.5: 47.926900000000003, 8.5: 47.870800000000003, 4: 48.2425, 6: 47.9298, 7: 47.983499999999999, 8: 47.896700000000003, 9: 47.833300000000001, 10: 47.8001, 11: 47.375100000000003, 12: 47.427199999999999, 13: 47.5471, 14: 47.6068, 15: 47.6616, 16: 47.679299999999998, 17: 47.735599999999998, 18: 48.030299999999997, 19: 48.091500000000003, 20: 48.142000000000003, 21: 48.188299999999998, 22: 48.271700000000003, 36.299999999999997: 47.075167934507682, 24: 48.3416, 11.199999999999999: 47.354999999999997, 27: 47.740299999999998, 28: 47.703400000000002, 29: 47.556800000000003, 30: 47.456499999999998, 31: 47.381100000000004, 32: 47.332900000000002, 33: 47.319800000000001, 34: 47.28636086132807, 35: 47.175599450858797, 36: 47.201012670573533, 35.600000000000001: 47.115630892576142, 11.1: 47.361764188029817, 403.10000000000002: 47.420621366868012, 35.5: 47.112494278018502, 'admiralty inlet': 48.173400000000001, 33.5: 47.323500000000003, 34.5: 47.195709921322582, 402.10000000000002: 47.376057851940473, 36.200000000000003: 47.146311826492713, 401: 47.490000000000002, 36.100000000000001: 47.167000000000002, 402: 47.363500000000002, 403: 47.40338787905921}
	
	lat = []
	lon = []
	for s in station_list:
		lat.append(lat_dict[s])
		lon.append(lon_dict[s])
	return np.array(lon), np.array(lat)

def high_res_station_to_lat_lon(station_list,n=5):
	lon, lat = station_to_lat_lon(station_list)
	
	hr_lat = []
	hr_lon = []
	
	for i in range(len(lat)-1):
		for j in range(n):
			frac = float(j)/float(n)
			hr_lat.append(lat[i]*(1-frac) + lat[i+1]*frac)
			hr_lon.append(lon[i]*(1-frac) + lon[i+1]*frac)
	hr_lat.append(lat[-1])
	hr_lon.append(lon[-1])
	return np.array(hr_lon), np.array(hr_lat)

def hood_canal_station_list():
#	return [24,22,21,20,19,18,7,7.5,8,8.5,9,10,17,16,15,14,13,401,12,11,402]
	return [7.5,8,8.5,9,10,17,16,15,14,13,401,12,11,11.1,11.2,402,402.1,403,403.1]
#	return [11,11.1,11.2,402,402.1,403,403.1] # This is a close up of the nasty bits at the end of Hood Canal that are difficult to get right

def main_basin_station_list():
	return [24,22,21,20,19,18,7,6,27,28,29,30,31,32,33,33.5,34,34.5,35,35.5,35.6,36, 36.1, 36.2, 36.3]
#	return [7,6,27,28,29,30,31,32,33,33.5,34]
	
def hood_canal_xy():
	return station_to_lat_lon(hood_canal_station_list())

def main_basin_xy():
	return station_to_lat_lon(main_basin_station_list())

def high_res_hood_canal_xy(n=1):
	if n == 1:
		return station_to_lat_lon(hood_canal_station_list())
	else:
		return high_res_station_to_lat_lon(hood_canal_station_list(),n)

def high_res_main_basin_xy(n=1):
	if n == 1:
		return station_to_lat_lon(main_basin_station_list())
	else:
		return high_res_station_to_lat_lon(main_basin_station_list(),n)
		
def latlon_to_km(lat1,lon1,lat2,lon2):
	RADIUS = 6378.0
#	d=2*asin(sqrt((sin((lat1-lat2)/2))^2 + cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))^2))
	
	lat1r = np.radians(lat1)
	lon1r = np.radians(lon1)
	lat2r = np.radians(lat2)
	lon2r = np.radians(lon2)	
	
	rads=2*(np.arcsin(np.sqrt(np.power(np.sin((lat1r-lat2r)/2.0),2.0) + np.cos(lat1r)*np.cos(lat2r)*np.power(np.sin((lon1r-lon2r)/2.0),2.0))))
	d = RADIUS*rads
	return d

def coords_to_km(coords):
	lon_list = coords['xm']
	lat_list = coords['ym']
	
	km_list = [0.0]
	for i in range(len(lon_list)-1):
		km_list.append(latlon_to_km(lat_list[i+1],lon_list[i+1],lat_list[i],lon_list[i]) + km_list[i])
	return km_list

def station_list_to_km(sl):
	lat,lon = station_to_lat_lon(sl)
	
	return coords_to_km({'xm':lon,'ym':lat})
	
def offset_calc(x1,y1,x2,y2,x3,y3,x4,y4):
	a = x2 - x1
	b = x3 - x4
	c = y2 - y1
	d = y3 - y4
	e = x3 - x1
	f = y3 - y1
	
	t = (d*e - b*f)/(a*d - b*c)
	
	if t >= 0.0 and t < 1.0:
		return t
	else:
		return None

def offset_region(coords,region='Admiralty Inlet'):
	if region == 'Admiralty Inlet':
		lon3 = -122.7090
		lat3 = 48.1931
		lon4 = -122.7774
		lat4 = 48.1267
	
	distances = coords_to_km(coords)
	lat = coords['ym']
	lon = coords['xm']
	
	rval = 0.0
	for i in range(len(lat)-1):
		lat1 = lat[i]
		lat2 = lat[i+1]
		lon1 = lon[i]
		lon2 = lon[i+1]
		t = offset_calc(lon1,lat1,lon2,lat2,lon3,lat3,lon4,lat4)
		if t >= 0.0 and t < 1.0:
			rval = distances[i] + t*(distances[i+1] - distances[i])
	return rval

def get_coastline(mode='regional'):
	module_dir = os.path.dirname(__file__)
	ncf_name = os.path.join(module_dir,'coastline.nc')
	try:
		ncf = nc.Dataset(ncf_name,mode='r')
		if mode.lower() == 'detailed':
			lat = ncf.variables['detailed_lat'][:]
			lon = ncf.variables['detailed_lon'][:]
			lat[lat==ncf.variables['detailed_lat']._FillValue] = np.nan
			lon[lon==ncf.variables['detailed_lon']._FillValue] = np.nan			
			
		else:
			lat = ncf.variables['regional_lat'][:]
			lon = ncf.variables['regional_lon'][:]
			lat[lat==ncf.variables['regional_lat']._FillValue] = np.nan
			lon[lon==ncf.variables['regional_lon']._FillValue] = np.nan
	except Exception, e:
		print(e)

	return lon,lat

def get_daves_section_var(section='AI_WB',var='label'):
	return daves_sections[section][var]

daves_numerical_vars = ['lon','lat','depth','depth_orig','along_dist']
daves_sections = {
	'AI_WB':{
		'lon':[-122.7270, -122.7122, -122.7029, -122.6899, -122.6770, -122.6677, -122.6510, -122.6417, -122.6324, -122.6306, -122.6269, -122.6250, -122.6232, -122.6213, -122.6158, -122.6102, -122.6046, -122.5972, -122.5898, -122.5805, -122.5638, -122.5490, -122.5305, -122.5119, -122.4971, -122.4841, -122.4748, -122.4711, -122.4637, -122.4637, -122.4674, -122.4656, -122.4544, -122.4396, -122.4192, -122.4025, -122.3858, -122.3729, -122.3654, -122.3580, -122.3506, -122.3413, -122.3339, -122.3284, -122.3209, -122.3209, -122.3265, -122.3376, -122.3432, -122.3488, -122.3580, -122.3766, -122.3914, -122.4081, -122.4229, -122.4378, -122.4526, -122.4656, -122.4860, -122.5082, -122.5230, -122.5305, -122.5397, -122.5453, -122.5509, -122.5564, -122.5564, -122.5546, -122.5546, -122.5527, -122.5466, -122.5359, -122.5252, -122.5123, -122.5029, -122.4934, -122.4896, -122.4913, -122.4990, -122.5067, -122.5059, -122.5020, -122.4943, -122.4861, -122.4780, -122.4686, -122.4591],
		'lat':[48.1723, 48.1636, 48.1549, 48.1449, 48.1337, 48.1237, 48.1124, 48.1025, 48.0925, 48.0787, 48.0650, 48.0525, 48.0413, 48.0313, 48.0176, 48.0051, 47.9939, 47.9839, 47.9727, 47.9627, 47.9577, 47.9515, 47.9465, 47.9402, 47.9327, 47.9277, 47.9215, 47.9140, 47.9015, 47.8903, 47.8766, 47.8641, 47.8541, 47.8516, 47.8566, 47.8616, 47.8678, 47.8753, 47.8841, 47.8978, 47.9090, 47.9203, 47.9340, 47.9465, 47.9614, 47.9739, 47.9864, 47.9989, 48.0139, 48.0288, 48.0426, 48.0488, 48.0563, 48.0663, 48.0762, 48.0850, 48.0950, 48.1025, 48.1074, 48.1149, 48.1237, 48.1374, 48.1511, 48.1649, 48.1786, 48.1936, 48.2085, 48.2223, 48.2322, 48.2435, 48.2525, 48.2626, 48.2687, 48.2762, 48.2834, 48.2923, 48.3033, 48.3134, 48.3215, 48.3270, 48.3350, 48.3405, 48.3443, 48.3486, 48.3521, 48.3555, 48.3578],
		'depth':[-72.1252, -67.9285, -72.5434, -78.7166, -81.8539, -84.8259, -117.4014, -128.7641, -144.3041, -154.0022, -130.6718, -125.2181, -120.7169, -116.4159, -115.7947, -96.9536, -90.2046, -91.9158, -104.0739, -107.2063, -101.3108, -106.9781, -115.7936, -121.8257, -146.3133, -131.3047, -130.5337, -141.6503, -156.0351, -170.0453, -194.9734, -214.7656, -217.6847, -194.3488, -169.1799, -154.1849, -160.7053, -180.6669, -204.5147, -205.0756, -189.0817, -178.2905, -176.9971, -180.1098, -175.7079, -160.4511, -158.5786, -145.7987, -132.3610, -138.7853, -130.4454, -138.8565, -128.4867, -118.2770, -110.7915, -107.6865, -105.5782, -118.4555, -105.5202, -107.4078, -96.8589, -86.2180, -85.3921, -83.2148, -83.8003, -81.2595, -74.8116, -71.3310, -69.2321, -64.2345, -56.4949, -41.8566, -27.1638, -18.5575, -15.1656, -14.6920, -13.9669, -10.0240, -8.0541, -7.2905, -4.6710, -4.3751, -4.6316, -4.7066, -4.7142, -4.7009, -4.6972],
		'depth_orig':[-80.2873, -68.2991, -75.5878, -77.0734, -98.6943, -70.3715, -130.5909, -128.2067, -156.3840, -171.6757, -137.8532, -132.0676, -127.8811, -119.2518, -126.2480, -98.5655, -90.9770, -96.4902, -108.5521, -111.9788, -102.0593, -114.6551, -126.5764, -125.4372, -157.2734, -126.9696, -137.6251, -147.4604, -155.5274, -171.4856, -203.3801, -225.7136, -252.2236, -202.7844, -177.0446, -149.9283, -160.4374, -178.6703, -214.4231, -206.6697, -189.7375, -179.2301, -178.1691, -184.1234, -175.9497, -160.2594, -159.8436, -144.3667, -159.5044, -153.1759, -166.1313, -170.7997, -150.8116, -120.0923, -112.9667, -108.0359, -111.5602, -123.2527, -118.7864, -112.4733, -96.9968, -86.4042, -86.1316, -84.6799, -86.3249, -83.4469, -76.3667, -73.1283, -72.2177, -75.0186, -65.4121, -64.1133, -29.0335, -21.7456, -18.0754, -17.7586, -22.0707, -8.8444, -3.8355, -1.3327, -0.0500, -0.0500, -0.0500, 'nan', -0.0500, 'nan', 'nan'],
		'along_dist':[0, 1.4666, 2.6559, 4.1244, 5.7005, 7.0057, 8.7635, 10.0689, 11.3743, 12.9058, 14.4558, 15.8493, 17.1049, 18.2228, 19.8032, 21.2502, 22.5650, 23.8039, 25.1684, 26.4748, 27.8348, 29.1384, 30.6257, 32.1703, 33.5530, 34.6672, 35.6456, 36.5223, 38.0150, 39.2630, 40.8132, 42.2068, 43.5919, 44.7323, 46.3511, 47.7133, 49.1376, 50.4136, 51.5306, 53.1530, 54.5179, 55.9442, 57.5664, 59.0136, 60.7668, 62.1535, 63.6006, 65.2153, 66.9300, 68.6446, 70.3183, 71.8604, 73.2411, 74.9044, 76.4677, 77.9357, 79.4986, 80.7715, 82.3834, 84.2320, 85.6994, 87.3209, 88.9940, 90.5741, 92.1542, 93.8685, 95.5325, 97.0640, 98.1734, 99.4289, 100.5289, 101.9034, 102.9438, 104.2091, 105.2719, 106.4864, 107.7381, 108.8677, 109.9316, 110.7662, 111.6665, 112.3392, 113.0455, 113.8159, 114.5302, 115.3261, 116.0685],
		'PRISM_ctd_ind':[3, 8, 13, 17, 26, 35, 55, 66],
		'label':['Adm. Inlet', 'Triple Jxn.', 'Saratoga Pass.', 'Skagit R.'],
		'label_ind':[3, 30, 55, 85]
	},
	
	'AI_HC':{
		'lon':[-122.7026, -122.6894, -122.6812, -122.6680, -122.6565, -122.6466, -122.6400, -122.6318, -122.6285, -122.6285, -122.6268, -122.6252, -122.6219, -122.6202, -122.6219, -122.6235, -122.6285, -122.6285, -122.6285, -122.6301, -122.6334, -122.6350, -122.6383, -122.6367, -122.6334, -122.6301, -122.6235, -122.6136, -122.6037, -122.5971, -122.6004, -122.6136, -122.6252, -122.6367, -122.6515, -122.6647, -122.6763, -122.6894, -122.7010, -122.7158, -122.7257, -122.7323, -122.7356, -122.7405, -122.7488, -122.7554, -122.7620, -122.7620, -122.7570, -122.7620, -122.7718, -122.7916, -122.8097, -122.8246, -122.8378, -122.8592, -122.8716, -122.8859, -122.9002, -122.9145, -122.9288, -122.9366, -122.9405, -122.9496, -122.9652, -122.9795, -122.9899, -122.9977, -123.0042, -123.0107, -123.0198, -123.0315, -123.0393, -123.0432, -123.0536, -123.0640, -123.0731, -123.0822, -123.0913, -123.0991, -123.1056, -123.1108, -123.1173, -123.1225, -123.1277, -123.1290, -123.1290, -123.1277, -123.1160, -123.1030, -123.0900, -123.0744, -123.0614],
		'lat':[48.1539, 48.1450, 48.1339, 48.1228, 48.1139, 48.1051, 48.0973, 48.0907, 48.0807, 48.0707, 48.0607, 48.0507, 48.0418, 48.0330, 48.0219, 48.0130, 48.0052, 47.9975, 47.9908, 47.9831, 47.9731, 47.9653, 47.9553, 47.9453, 47.9331, 47.9209, 47.9110, 47.9054, 47.8988, 47.8866, 47.8755, 47.8655, 47.8577, 47.8500, 47.8411, 47.8333, 47.8244, 47.8167, 47.8100, 47.8012, 47.7912, 47.7790, 47.7679, 47.7546, 47.7446, 47.7346, 47.7246, 47.7124, 47.7013, 47.6913, 47.6836, 47.6780, 47.6736, 47.6703, 47.6680, 47.6625, 47.6542, 47.6463, 47.6393, 47.6323, 47.6253, 47.6165, 47.6078, 47.6008, 47.5964, 47.5833, 47.5728, 47.5623, 47.5544, 47.5474, 47.5369, 47.5255, 47.5142, 47.5028, 47.4923, 47.4809, 47.4722, 47.4625, 47.4547, 47.4459, 47.4363, 47.4267, 47.4162, 47.4066, 47.3969, 47.3873, 47.3777, 47.3716, 47.3689, 47.3672, 47.3672, 47.3663, 47.3619],
		'depth':[-73.1053, -77.6711, -90.7865, -87.0021, -113.7509, -128.2548, -145.2415, -147.4307, -150.9238, -137.9719, -127.4332, -124.6678, -119.9143, -115.2329, -117.5180, -109.2843, -109.0942, -103.6687, -92.9045, -77.2303, -66.7233, -65.4203, -76.9665, -91.6236, -93.4275, -89.0655, -90.3293, -91.7785, -93.8039, -93.0508, -86.6447, -74.4546, -68.7260, -63.9077, -59.4231, -53.8171, -54.5158, -53.5801, -52.5939, -49.9461, -62.0773, -73.0567, -76.8136, -85.8187, -92.5836, -93.3260, -88.3262, -81.9592, -76.1031, -83.5185, -88.2079, -96.9549, -95.7813, -99.5001, -116.8061, -129.3801, -133.9035, -144.5435, -138.1318, -139.5590, -140.7119, -152.7918, -159.6839, -153.1871, -142.6536, -141.3210, -135.5897, -141.6485, -136.3260, -124.0167, -84.7916, -99.6415, -99.9918, -97.2696, -95.7214, -88.0845, -93.1406, -92.8867, -89.0915, -90.6100, -90.7535, -88.0933, -86.4537, -81.8304, -69.7217, -66.2787, -69.5117, -66.4397, -42.7274, -25.3152, -23.7716, -20.7440, -14.9898],
		'depth_orig':[-74.6230, -74.4055, -110.2748, -75.7537, -130.2792, -135.9251, -148.4656, -161.4767, -170.4817, -147.6731, -133.5608, -129.9250, -125.3487, -120.1464, -119.8713, -105.7969, -110.7673, -107.5087, -94.9731, -78.2472, -67.1327, -61.8776, -78.9527, -105.6035, -114.1106, -96.4530, -95.3211, -100.4642, -115.5858, -114.0148, -95.6119, -84.9722, -81.8498, -73.4327, -65.4853, -60.8050, -61.0048, -56.5926, -54.6884, -52.5331, -69.8933, -77.3357, -84.6057, -96.0947, -115.0989, -108.0324, -95.4909, -87.2057, -77.7738, -95.7063, -103.3177, -113.5022, -100.0551, -99.3995, -125.1088, -128.9849, -135.0734, -150.0098, -156.2140, -152.7505, -155.4808, -166.0193, -168.7131, -163.0757, -158.1357, -162.0572, -158.2650, -157.8033, -148.4720, -144.0830, -100.9106, -131.5868, -136.0433, -137.6791, -135.0406, -129.6500, -144.8182, -119.7838, -118.1513, -121.6619, -121.5995, -116.3518, -112.2437, -112.9469, -94.4990, -91.8128, -87.1988, -76.0735, -55.0406, -39.6097, -34.1274, -25.4283, -15.6809],
		'along_dist':[0, 1.3884, 2.7641, 4.3375, 5.6432, 6.8722, 7.8640, 8.8237, 9.9597, 11.0690, 12.1850, 13.3010, 14.3170, 15.3106, 16.5492, 17.5428, 18.4806, 19.3434, 20.0829, 20.9544, 22.0905, 22.9619, 24.0980, 25.2140, 26.5919, 27.9697, 29.1828, 30.1432, 31.1870, 32.6290, 33.8858, 35.3678, 36.5861, 37.8045, 39.2863, 40.5946, 41.9034, 43.2119, 44.3468, 45.8293, 47.1618, 48.6041, 49.8610, 51.3855, 52.6542, 53.8680, 55.0818, 56.4376, 57.7244, 58.8937, 60.0301, 61.6329, 63.0762, 64.2460, 65.2630, 66.9807, 68.2920, 69.6746, 70.9980, 72.3214, 73.6450, 74.7793, 75.7945, 76.8288, 78.0948, 79.9047, 81.3079, 82.6129, 83.6146, 84.5325, 85.8842, 87.4231, 88.8160, 90.1134, 91.5172, 93.0030, 94.1914, 95.4607, 96.5712, 97.7064, 98.8822, 100.0209, 101.2858, 102.4245, 103.5633, 104.6372, 105.7067, 106.3942, 107.3218, 108.3194, 109.2978, 110.4760, 111.5686],
		'PRISM_ctd_ind':[1, 7, 13, 19, 28, 35, 39, 45, 49, 55, 62, 69, 74, 81, 87],
		'label':['Adm. Inlet', 'Sill', 'Dabob Bay', 'Hoodsport', 'Lynch Cove'],
		'label_ind':[3, 35, 55, 81, 92]
	},
	
	'JdF_SoG':{
		'lon':[-123.0193, -123.0042, -122.9992, -122.9992, -123.0067, -123.0143, -123.0218, -123.0420, -123.0496, -123.0622, -123.0723, -123.0874, -123.1025, -123.1202, -123.1353, -123.1504, -123.1655, -123.1731, -123.1807, -123.1857, -123.1908, -123.1983, -123.2084, -123.2160, -123.2185, -123.2185, -123.2210, -123.2210, -123.2235, -123.2311, -123.2412, -123.2538, -123.2563, -123.2513, -123.2437, -123.2261, -123.2084, -123.1908, -123.1706, -123.1529, -123.1353, -123.1151, -123.0975, -123.0798, -123.0597, -123.0420, -123.0319, -123.0218, -123.0092, -123.0042, -123.0092, -123.0118, -123.0319, -123.0571, -123.0874, -123.1151, -123.1403, -123.1655, -123.1857, -123.2059, -123.2210, -123.2387, -123.2563, -123.2765, -123.2916, -123.3092, -123.3244, -123.3395, -123.3546, -123.3647, -123.3748, -123.3874, -123.3748, -123.3521, -123.3294, -123.3017, -123.2714, -123.2462, -123.2210, -123.2008, -123.1731, -123.1529],
		'lat':[48.3339, 48.3473, 48.3639, 48.3805, 48.3938, 48.4055, 48.4155, 48.4288, 48.4388, 48.4488, 48.4554, 48.4637, 48.4704, 48.4754, 48.4820, 48.4920, 48.5020, 48.5136, 48.5253, 48.5353, 48.5469, 48.5569, 48.5686, 48.5785, 48.5902, 48.6052, 48.6168, 48.6285, 48.6401, 48.6534, 48.6634, 48.6767, 48.6900, 48.7000, 48.7083, 48.7150, 48.7183, 48.7183, 48.7200, 48.7266, 48.7300, 48.7349, 48.7399, 48.7449, 48.7532, 48.7599, 48.7715, 48.7832, 48.7948, 48.8065, 48.8198, 48.8348, 48.8414, 48.8448, 48.8464, 48.8514, 48.8547, 48.8614, 48.8697, 48.8797, 48.8913, 48.9030, 48.9146, 48.9263, 48.9396, 48.9529, 48.9662, 48.9795, 48.9928, 49.0078, 49.0211, 49.0361, 49.0561, 49.0694, 49.0794, 49.0894, 49.0960, 49.1027, 49.1110, 49.1160, 49.1193, 49.1226],
		'depth':[-145.6121, -154.9663, -154.9064, -152.1677, -155.7285, -162.7622, -172.1557, -180.6019, -177.1227, -180.6902, -188.5748, -208.3572, -194.1927, -184.8497, -229.7557, -270.7973, -258.7998, -244.0210, -243.2408, -240.9671, -242.9330, -251.7905, -262.9741, -264.1829, -256.8572, -248.8137, -244.5476, -242.9752, -239.8129, -239.1362, -258.5232, -295.2163, -310.2610, -305.4461, -273.2152, -208.7353, -167.5729, -179.3427, -208.1616, -191.6202, -155.9835, -134.9392, -133.9871, -138.4068, -174.1352, -181.9742, -169.5464, -152.9165, -174.1030, -188.8982, -199.5301, -200.6982, -195.5081, -193.5516, -185.5922, -172.2676, -161.7746, -156.6017, -156.0640, -152.9986, -155.1711, -155.3900, -160.1532, -168.3474, -175.9048, -186.6216, -197.5065, -210.0864, -227.1121, -241.4641, -251.7601, -261.5386, -238.9857, -197.4275, -140.2425, -42.0877, -7.8559, -5.3316, -6.2845, -8.2251, -8.0384, -7.8511],
		'depth_orig':[-141.9237, -154.5394, -161.3531, -153.9435, -155.3613, -166.6722, -175.8714, -187.5720, -188.0398, -195.6347, -184.1267, -221.8005, -187.0240, -131.2245, -222.0819, -285.3259, -270.3345, -218.7556, -248.1045, -247.8425, -246.0096, -260.7029, -269.9102, -267.0873, -263.9724, -253.3104, -247.6621, -243.4448, -246.9784, -242.1123, -290.6214, -312.4277, -330.0115, -355.9161, -324.9971, -225.8884, -215.6733, -256.6477, -247.9950, -228.0047, -168.5173, -134.5303, -137.4053, -121.4901, -187.9220, -205.1471, -193.3375, -173.5089, -182.0094, -168.4780, -219.5183, -192.8690, -193.6182, -196.5662, -201.0859, -168.9836, -154.9015, -153.8517, -151.9308, -153.3023, -148.2517, -152.7498, -154.1578, -165.5843, -172.5137, -185.9609, -197.1376, -208.6338, -228.2415, -243.0968, -252.3216, -264.8352, -244.4368, -199.9745, -138.9301, -22.3451, -7.0008, -5.3829, -3.9514, -2.4842, -0.7844, 'nan'],
		'along_dist':[0, 1.8537, 3.7397, 5.5886, 7.1695, 8.5789, 9.8207, 11.9182, 13.1598, 14.6069, 15.6554, 17.1035, 18.4411, 19.8545, 21.1918, 22.7640, 24.3360, 25.7450, 27.1538, 28.3236, 29.6699, 30.9110, 32.4026, 33.6435, 34.9510, 36.6150, 37.9224, 39.2167, 40.5241, 42.1040, 43.4376, 45.1822, 46.6728, 47.8422, 48.9203, 50.4107, 51.7563, 53.0501, 54.5401, 56.0303, 57.3756, 58.9544, 60.3616, 61.7687, 63.5118, 65.0012, 66.4913, 67.9814, 69.5710, 70.9168, 72.4412, 74.1155, 75.7655, 77.6459, 79.8659, 81.9682, 83.8482, 85.8342, 87.5743, 89.4191, 91.1211, 92.9478, 94.7743, 96.7349, 98.5807, 100.5420, 102.3875, 104.2328, 106.0779, 107.8971, 109.5487, 111.4494, 113.8505, 116.0680, 118.0576, 120.3606, 122.6828, 124.6605, 126.7143, 128.2826, 130.3331, 131.8457],
		'IOS_ctd_ind':[15, 22, 27, 35, 40, 46, 57, 72],
		'JEMS_ctd_ind':[0, 6],
		'label':['Juan de Fuca', 'Haro St.', 'Bndy. Pass', 'Fraser R.'],
		'label_ind':[3, 22, 40, 79]
	}
}


for key in daves_sections.keys():
	for var in daves_numerical_vars:
		daves_sections[key][var] = np.array(daves_sections[key][var],'double')