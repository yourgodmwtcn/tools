import os

import netCDF4 as nc
import numpy as np

import load_grid
#import plot_utils
import utils

__version__ = '0.1.6'

def extract_from_file(file='',varname='zeta',extraction_type='full',**kwargs):
	# Check the netCDF file for existence and if the variable is in it
	if not os.path.exists(file):
		raise IOError('File %s could not be located on the filesystem' %file)
	ncf = nc.Dataset(file,mode='r')
	
	if varname not in ncf.variables:
		raise IOError('File %s does not have a variable named %s' % (file, varname))
	
	# start getting data
	ncvar = ncf.variables[varname]
	dims = ncvar.dimensions
	ndims = len(dims)
	shape = ncvar.shape
#	print('var: %s, dims: %s, shape: %s' %(varname, str(dims), str(shape)))
	if not ndims == 3 and not ndims == 4:
		raise TypeError('ndims is neither 3 nor 4')
	if not dims[0] == 'ocean_time':
		raise TypeError('first dimension is not ocean_time')
	if not shape[0] == 1:
		raise TypeError('first dimension is not of length one')
	
	grid = load_grid.load_grid(file)
	coords = {}
	if ndims == 3:
		if dims[1] == 'eta_rho' and dims[2] == 'xi_rho':
#			print('G.lat %s' % str(grid['lat']))
#			print('G.lon %s' % str(grid['lon']))
			y2 = grid['lat'][:]
			x2 = grid['lon'][:]
			mask2 = grid['mask'][:]
		elif dims[1] == 'eta_u' and dims[2] == 'xi_u':
#			print('G.latu %s' % str(grid['lat_u']))
#			print('G.lonv %s' % str(grid['lon_u']))
			y2 = grid['latu'][:]
			x2 = grid['lonu'][:]
			mask2 = grid['masku'][:]
		elif dims[1] == 'eta_v' and dims[2] == 'xi_v':
#			print('G.latv %s' % str(grid['lat_v']))
#			print('G.lonv %s' % str(grid['lon_v']))
			y2 = grid['latv'][:]
			x2 = grid['lonv'][:]
			mask2 = grid['maskv'][:]
		else:
			raise TypeError('Unable to determine which gird to use')
		
		if extraction_type == 'full' or extraction_type == 'surface':
#			data = np.squeeze(ncvar[:])
#			data[mask2==0] = np.NaN
			data = np.ma.array(np.squeeze(ncvar[:]),mask=(mask2==0))
			
			coords['ym'] = y2
			coords['xm'] = x2
		
		if (extraction_type == 'profile' or extraction_type == 'profiles' or
		    extraction_type == 'point' or extraction_type == 'points') and \
		   (kwargs.has_key('y') and kwargs.has_key('x')):
			print('In profiles')
			xm = np.array(kwargs['x'])
			ym = np.array(kwargs['y'])
			if not (xm.ndim == ym.ndim) or not (xm.shape == ym.shape):
				if xm.size == 1:
					xm = np.repmat(xm,ym.shape)
				elif ym.size == 1:
					ym = np.repmat(ym,xm.shape)
#				else:
#					raise RuntimeError('The x and y chosen to extract a point or profile on this 2D variable are incompatible.')

			data2 = np.ma.array(np.squeeze(ncvar[:]),mask=(mask2==0))
#			data = utils.interp_2d(lat=y2,lon=x2,data=data2,lati=ym,loni=xm)
			mask = utils.interp_2d_xy(y=y2,x=x2,data=mask2,yi=ym,xi=xm)
			data = utils.interp_2d_xy(y=y2,x=x2,data=data2,yi=ym,xi=xm)
			data = np.ma.array(data,mask=(mask < 1))
			coords['ym'] = ym
			coords['xm'] = xm

	if ndims == 4:
		xm = []
		ym = []
		zm = []
		
		K = shape[1]
		J = shape[2]
		I = shape[3]
		
		lat = grid['lat'][:]
		lon = grid['lon'][:]
		
		if dims[2] == 'eta_rho' and dims[3] == 'xi_rho':
			y2 = lat
			x2 = lon
			mask2 = grid['mask'][:]
		
		elif dims[2] == 'eta_u' and dims[3] == 'xi_u':
			y2 = grid['latu'][:]
			x2 = grid['lonu'][:]
			mask2 = grid['masku'][:]
		
		elif dims[2] == 'eta_v' and dims[3] == 'xi_v':
			y2 = grid['latv'][:]
			x2 = grid['lonv'][:]
			mask2 = grid['maskv'][:]
		else:
			raise TypeError('Unable to determine which gird to use. dims[2] = %s, dims[3] = %s' % (dims[2], dims[3]))
		
		if dims[1] == 's_rho':
			cs = grid['cs'][:]
		elif dims[1] == 's_w':
			cs = grid['csw'][:]
		elif K==1:
			cs = 0
		else:
			raise TypeError('Unable to determine which cs to use. dim[1] = %s' % dim[1])
		
		data = []
		if extraction_type == 'full':
			# get zeta		
			try:
				if dims[2] == 'eta_rho' and dims[3] == 'xi_rho':
					zeta2 = ncf.variables['zeta'][:]
				else:
					zeta2 = utils.interp_2d_xy(y=lat,x=lon,data=ncf.variables['zeta'][:],yi=y2,xi=x2)
			except Exception, e:
				print(e)
				zeta2 = np.zeros((len(y2),len(x2)))
			zeta2[zeta2>1000] = 0
	
			# get H
			if dims[2] == 'eta_rho' and dims[3] == 'xi_rho':
				H2 = grid['H'][:]
			else:
				H2 = utils.interp_2d_xy(y=lat,x=lon,data=grid['H'][:],yi=y2,xi=x2)
			
			x3 = np.tile(x2.reshape(1, J, I),(K,1,1))
			y3 = np.tile(y2.reshape(1, J, I),(K,1,1))
			mask3 = np.tile(mask2.reshape(1, J, I),(K,1,1))
			zeta3 = np.tile(zeta2.reshape(1, J, I),(K,1,1))
			H3 = np.tile(H2.reshape(1, J, I),(K,1,1))
			cs3 = np.tile(cs.reshape(K,1,1),(1,J,I))
			z3 = zeta3 + cs3*(zeta3 + H3)
			
			zm = z3
			ym = y3
			xm = x3
			data = np.ma.array(np.squeeze(ncvar[:]),mask=(mask3==0))
		
		elif extraction_type == 'surface':
			# get zeta		
			try:
				if dims[2] == 'eta_rho' and dims[3] == 'xi_rho':
					zeta2 = ncf.variables['zeta'][:]
				else:
					zeta2 = utils.interp_2d_xy(y=lat,x=lon,data=ncf.variables['zeta'][:],yi=y2,xi=x2)
			except Exception, e:
				print(e)
				zeta2 = np.zeros((len(y2),len(x2)))
			zeta2[zeta2>1000] = 0
			data = np.ma.array(np.squeeze(ncvar[0,K-1,:,:]),mask=(mask2==0))
			zm = zeta2
			ym = y2
			xm = x2
		
		elif (extraction_type == 'profile' or extraction_type == 'profiles') \
							and kwargs.has_key('x') and kwargs.has_key('y'):
			
			if dims[1] == 's_rho':
				K = K + 2
				cs_list = [-1.0]
				tmp = grid['cs'][:]
				for t in tmp:
					cs_list.append(t)
				cs_list.append(0.0)
				cs = np.array(cs_list)
			elif dims[1] == 's_w':
				cs = grid['csw'][:]
			elif K==1:
				cs = 0
			else:
				raise TypeError('Unable to determine which cs to use. dim[1] = %s' % dim[1])
			
			x3 = np.tile(x2.reshape(1, J, I),(K,1,1))
			y3 = np.tile(y2.reshape(1, J, I),(K,1,1))
			mask3 = np.tile(mask2.reshape(1, J, I),(K,1,1))
			cs3 = np.tile(cs.reshape(K,1,1),(1,J,I))
			
			# TODO: add a bunch of code that takes the x an y and makes the appropriate xi, yi, zi 3d matricies for interpolation with utils.interp_3d()
			x = kwargs['x']
			y = kwargs['y']
			
			if not x.shape == y.shape or not x.ndim == 1:
				raise(TypeError,'In profiles x and y must be vectors that are the same length.')

			if kwargs.has_key('cs'):
				csi = kwargs['cs']
				zi = np.tile(csi.reshape(csi.size,1),(1,x.size))
			else:
				zi = np.tile(cs.reshape(cs.size,1),(1,x.size))

			xi = np.tile(x.reshape(1,x.size),(zi.shape[0],1))
			yi = np.tile(y.reshape(1,y.size),(zi.shape[0],1))

			zetai = np.zeros(x.shape)
			Hi = np.zeros(x.shape)
			for i in range(len(x)):
				zetai[i] = utils.interp_2d(x=lon,y=lat,data=np.squeeze(ncf.variables['zeta'][:]),yi=y[i],xi=x[i])
				Hi[i] = utils.interp_2d(x=lon,y=lat,data=np.squeeze(grid['H'][:]),yi=y[i],xi=x[i])
			z = zetai + zi*(zetai + Hi)
			
			if dims[1] == 's_rho':
				rd = np.squeeze(ncvar[:])
				rd = np.concatenate((rd[0:1,:,:],rd,rd[-1:,:,:]))
			elif dims[1] == 's_w':
				rd = np.squeeze(ncvar[:])
			
			data = utils.interp_3d(x=x3,y=y3,z=cs3,data=rd,xi=xi,yi=yi,zi=zi)
#			data = np.ma.array(data,mask=(data>100))
			data = np.ma.array(data,mask=(z>100))
#			data[data>100] = 0.
			z[z>100] = 0.0
#			print(data)
			xm = x
			ym = y
			zm = z
			
		coords = {}
		if len(xm) > 0:
			coords['xm'] = xm
		if len(ym) > 0:
			coords['ym'] = ym
		if len(zm) > 0:
			coords['zm'] = zm
	grid['ncfile'].close()
	ncf.close()
	return (data,coords)