import numpy as np
from matplotlib.mlab import griddata

__version__ = '0.1'

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
	lat_dict = {}
	lon_dict = {}
	
	lat_dict[4] = 48.2425
	lon_dict[4] = -122.542
	
#	lat_dict[6] = 47.925
#	lon_dict[6] = -122.4685
	
	lat_dict[6] = 47.9298
	lon_dict[6] = -122.4932
	
	lat_dict[7] = 47.9835
	lon_dict[7] = -122.6201
	
	lat_dict[7.5] = 47.9269
	lon_dict[7.5] = -122.6418
	
	lat_dict[8] = 47.8967
	lon_dict[8] = -122.6053
	
	lat_dict[8.5] = 47.8708
	lon_dict[8.5] = -122.5848
	
	lat_dict[9] = 47.8333
	lon_dict[9] = -122.6673
	
	lat_dict[10] = 47.8001
	lon_dict[10] = -122.7198

	lat_dict[11] = 47.3751
	lon_dict[11] = -123.1375
	
	lat_dict[11.1] = 47.36176418802982
	lon_dict[11.1] = -123.063768617271
	
	lat_dict[11.2] = 47.3550
	lon_dict[11.2] = -123.0305
	
	lat_dict[12] = 47.4272
	lon_dict[12] = -123.1142

	lat_dict[13] = 47.5471
	lon_dict[13] = -123.008

	lat_dict[14] = 47.6068
	lon_dict[14] = -122.9399

	lat_dict[15] = 47.6616
	lon_dict[15] = -122.8601
	
# true station 16 lat lon
#	lat_dict[16] = 47.6917
#	lon_dict[16] = -122.6074

# fake station 16 lat lon
	lat_dict[16] = 47.6793
	lon_dict[16] = -122.7578
	
	lat_dict[17] = 47.7356
	lon_dict[17] = -122.7614

	lat_dict[18] = 48.0303
	lon_dict[18] = -122.6169

	lat_dict[19] = 48.0915
	lon_dict[19] = -122.6318

	lat_dict[20] = 48.142
	lon_dict[20] = -122.6848

	lat_dict[21] = 48.1883
	lon_dict[21] = -122.8504

	lat_dict[22] = 48.2717
	lon_dict[22] = -123.0189

	lat_dict[24] = 48.3416
	lon_dict[24] = -123.1249

#	lat_dict[27] = 47.8133
#	lon_dict[27] = -122.4050
	
	lat_dict[27] = 47.7403
	lon_dict[27] = -122.4103
	
	lat_dict[28] = 47.7034
	lon_dict[28] = -122.4544

	lat_dict[29] = 47.5568
	lon_dict[29] = -122.4433

	lat_dict[30] = 47.4565
	lon_dict[30] = -122.4084

#	lat_dict[31] = 47.3937
#	lon_dict[31] = -122.3601

	lat_dict[31] = 47.3811
	lon_dict[31] = -122.3487

	lat_dict[32] = 47.3329
	lon_dict[32] = -122.4427

	lat_dict[33] = 47.3198
	lon_dict[33] = -122.5008
	
	lat_dict[33.5] = 47.3235
	lon_dict[33.5] = -122.5621
	
	lat_dict[34] = 47.28636086132807
	lon_dict[34] = -122.5372204355572
	
	lat_dict[34.5] = 47.19570992132258
	lon_dict[34.5] = -122.6043449761104
	
	lat_dict[35] = 47.1755994508588
	lon_dict[35] = -122.6519213323306
	
	lat_dict[35.5] = 47.1124942780185
	lon_dict[35.5] = -122.6858130539733
	
	lat_dict[35.6] = 47.11563089257614
	lon_dict[35.6] = -122.7313463437406
	
	lat_dict[36] = 47.20101267057353
	lon_dict[36] = -122.825703220149
	
	lat_dict[36.1] = 47.1670
	lon_dict[36.1] = -122.8573269666242
	
	lat_dict[36.2] = 47.14631182649271
	lon_dict[36.2] = -122.9157302320745
	
	lat_dict[36.3] = 47.07516793450768
	lon_dict[36.3] = -122.9127057495704
	
	lat_dict[401] = 47.49
	lon_dict[401] = -123.0567

	lat_dict[402] = 47.3635
	lon_dict[402] = -123.0167366137608
	
	lat_dict[402.1] = 47.37605785194047
	lon_dict[402.1] = -123.0013100422375
	
	lat_dict[403] = 47.40338787905921
	lon_dict[403] = -122.9112721092632
	
	lat_dict[403.1] = 47.42062136686801
	lon_dict[403.1] = -122.8790276962233
	
	lat_dict['admiralty inlet'] = 48.1734
	lon_dict['admiralty inlet'] = -122.7557
	
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
			