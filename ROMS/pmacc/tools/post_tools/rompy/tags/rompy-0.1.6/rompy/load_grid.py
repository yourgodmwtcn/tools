# module load_grid
import netCDF4 as nc

import numpy as np
#from scipy.interpolate.fitpack2 import  SmoothBivariateSpline

__version__ = '0.1.6'

def load_grid(files):
	'''
load_grid returns a dictionary of netCDF4 variables. The actual extraction of the data in them from
the netCDF file must be done by the user, hopefully at the last minute, so as to reduce the amount
of wasted memory.
	'''
	if files.__class__ == str:
		f = files
	elif files.__class__ == list:
		f = files[0]
	else:
		raise TypeError('wrong type sent to load_grid')
		return False
	
	ncf = nc.Dataset(f,mode='r')
	
	G = {}
	G['ncfile'] = ncf
	G['lon'] = ncf.variables['lon_rho']
	G['lat'] = ncf.variables['lat_rho']
	G['lonu'] = ncf.variables['lon_u']
	G['latu'] = ncf.variables['lat_u']
	G['lonv'] = ncf.variables['lon_v']
	G['latv'] = ncf.variables['lat_v']
	G['cs'] = ncf.variables['Cs_r']
	G['csw'] = ncf.variables['Cs_w']
	
	G['mask'] = ncf.variables['mask_rho']
	G['masku'] = ncf.variables['mask_u']
	G['maskv'] = ncf.variables['mask_v']
	
	G['H'] = ncf.variables['h']
#	G['Hu'] = griddata(G['lat'].reshape(G['lat'].size),G['lon'].reshape(G['lon'].size),G['H'].reshape(G['H'].size),G['latu'],G['lonu'])
#	G['Hv'] = griddata(G['lat'].reshape(G['lat'].size),G['lon'].reshape(G['lon'].size),G['H'].reshape(G['H'].size),G['latv'],G['lonv'])
	return G
