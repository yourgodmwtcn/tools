#!/usr/bin/env python
import os
from optparse import OptionParser
from subprocess import Popen,PIPE

import netCDF4 as nc
import numpy as np

def dim2vals(dims,ncf):
	dim_vals = []
	for i in range(len(dims)):
		dim_vals.append(len(ncf.dimensions[dims[i]]))
	return tuple(dim_vals)

def z_coords(N=20, theta_s=4.0, theta_b=0.5,tcline=0.0):
	# it is important to note that N here is to be the same size as s_rho and Cs_r. S_w and Cs_w will be of length N+1
	if not theta_s == 0:
		cff1 = 1.0/np.sinh(theta_s)
		cff2 = 0.5/np.tanh(0.5*theta_s)
	
	s_rho = (np.linspace(-N+1,0,N) - 0.5)/N
	s_w = (np.linspace(-N,0,N+1))/N
	
	if not theta_s == 0:
		Cs_r = (1.0-theta_b)*cff1*np.sinh(theta_s*s_rho) + theta_b*( cff2*np.tanh(theta_s*(s_rho + 0.5)) - 0.5 )
		Cs_w = (1.0-theta_b)*cff1*np.sinh(theta_s*s_w) + theta_b*( cff2*np.tanh(theta_s*(s_w + 0.5)) - 0.5 )
	
	S = {'N':N, 's_rho':s_rho, 's_w':s_w, 'Cs_r':Cs_r, 'Cs_w':Cs_w}
	return S

def new_var_copy_of_old_var(ncf, ovn, nvn, nc_var_type='f8'):
	ov = ncf.variables[ovn]
	ovd = ov.dimensions
	
	ovdv = dim2vals(ovd, ncf)
	
	try:
		nv = ncf.createVariable(nvn, nc_var_type, ovd)
	except RuntimeError, e:
		print(e)
		print('Overwriting variable "%s"' % nvn)
		nv = ncf.variables[nvn]
	
	n = nv[:]
	o = ov[:]
	for i in range(n.size):
		n.flat[i] = o.flat[i]
	nv[:] = n[:]

def uniform_new_var_from_old(ncf, ovn, nvn, fill_val=0.0, long_name=None, units=None,nc_var_type='f8'):
	ov = ncf.variables[ovn]
	ovd = ov.dimensions
	
	ovdv = dim2vals(ovd, ncf)
	
	try:
		nv = ncf.createVariable(nvn, nc_var_type, ovd)
	except RuntimeError, e:
		print(e)
		print('Overwriting variable "%s"' % nvn)
		nv = ncf.variables[nvn]
	
	nv[:] = fill_val*np.ones(ovdv)
	
#	copy all attributes from the old variable to the new one. If that is a bad idea use the commented code above to hard code which ones to copy.
	for attr in ov.ncattrs():
		nv.__setattribute__(attr)
	
	if long_name is not None:
		nv.long_name = long_name
	
	if units is not None:
		nv.units = units
	
#	for attr in nv.ncattrs():
#		print('%s: %s' % (attr, nv.__getattribute__(attr)))

def new_var_linear_function_of_old_var(ncf, ovn, nvn, a=0.0, b=1.0, long_name=None, units=None,nc_var_type='f8'):
#	the new variable is a linear function of the old or base variable.
#	nv = b*ov + a.
#	Default behaviour is to copy the old variable to the new variable.
	ov = ncf.variables[ovn]
	ovd = ov.dimensions
	
	ovdv = dim2vals(ovd, ncf)
	
	try:
		nv = ncf.createVariable(nvn, nc_var_type, ovd)
	except RuntimeError, e:
		print(e)
		print('Overwriting variable "%s"' % nvn)
		nv = ncf.variables[nvn]
	
	n = nv[:]
	o = ov[:]
	for i in range(n.size):
		n.flat[i] = a + b*o.flat[i]
	nv[:] = n[:]
	if long_name is not None:
		nv.long_name = long_name
	
	if units is not None:
		nv.units = units

def new_var_function_of_old_var_and_depth(ncf,ovn,nvn,grd_ncf, long_name=None, units=None,nc_var_type='f8'):
	h_var = grd_ncf.variables['h']
	
	suffix = ovn.split('_')[-1]
	
	if suffix == 'east':
		h = h_var[:,-1]
	elif suffix == 'west':
		h = h_var[:,0]
	elif suffix == 'north':
		h = h_var[-1,:]
	elif suffix == 'south':
		h = h_var[0,:]
	else:
		h = h_var[:]
	
	S = z_coords(N=len(ncf.dimensions['s_rho']))
	
	if h.ndim == 1:
		z = np.outer(S['Cs_r'],h)
		print(h.shape)
		print(S['Cs_r'].shape)
		print(z.shape)
		
	elif h.ndim == 2:
		z = np.zeros((S['N'],h.shape[0],h.shape[1]))
		for i in range(S['N']):
			z[i,:,:] = S['Cs_r'][i]*h

	ov = ncf.variables[ovn]
	ovd = ov.dimensions
	
	ovdv = dim2vals(ovd, ncf)
	print(ovdv)
	
	try:
		nv = ncf.createVariable(nvn, nc_var_type, ovd)
	except RuntimeError, e:
		print(e)
		print('Overwriting variable "%s"' % nvn)
		nv = ncf.variables[nvn]
	
	n = nv[:]
	o = ov[:]
	
	if suffix == 'east':
		salinity = ncf.variables['salt_east'][:]
	elif suffix == 'west':
		salinity = ncf.variables['salt_west'][:]
	elif suffix == 'north':
		salinity = ncf.variables['salt_north'][:]
	elif suffix == 'south':
		salinity = ncf.variables['salt_south'][:]
	else:
		salinity = ncf.variables['salt'][:]
	
	if len(ovdv) == 4:
		# assume var has dims like (time,z,lat,lon)
		for t in range(ovdv[0]):
			for k in range(ovdv[1]):
				for j in range(ovdv[2]):
					for i in range(ovdv[3]):
						if np.abs(z[k,0,i]) < 0.001:
							n[t,k,j,i] = z[k,j,i]*salinity[t,k,j,i]
						else:
							n[t,k,j,i] = z[k,j,i]*salinity[t,k,j,i]/z[k,0,i]
	elif len(ovdv) == 3:
		# assume var has dims like (time, z, lat) or (time,z,lon)
		for t in range(ovdv[0]):
			for j in range(ovdv[1]):
				for i in range(ovdv[2]):
					if np.abs(z[0,i]) < 0.001:
						n[t,j,i] = z[j,i]*salinity[t,j,i]
					else:
						n[t,j,i] = z[j,i]*salinity[t,j,i]/z[0,i]
	
#	for i in range(n.size):
#		n.flat[i] = a + b*o.flat[i]
	nv[:] = n[:]
	
	
	if long_name is not None:
		nv.long_name = long_name
	
	if units is not None:
		nv.units = units

if __name__ == '__main__':
	parser = OptionParser()
	
	parser.add_option('-t','--test',
						dest='testing',
						action='store_true',
						default=True,
						help='Testing mode will copy t he file first and edit the copy only.')
	
	(options, args) = parser.parse_args()
	
	bry_file_name = 'dataset/ocean_bry.nc'
	ini_file_name = 'dataset/ocean_ini.nc'
	clm_file_name = 'dataset/ocean_clm.nc'
	grid_file_name = 'dataset/grid.nc'
	
	if options.testing:
		p = Popen('cp %s test_bry.nc' % bry_file_name, shell=True, stdout=PIPE)
		sts = os.waitpid(p.pid, 0)[1]
		bry_file_name = 'test_bry.nc'
		
		p = Popen('cp %s test_ini.nc' % ini_file_name, shell=True, stdout=PIPE)
		sts = os.waitpid(p.pid, 0)[1]
		ini_file_name = 'test_ini.nc'
		
#		p = Popen('cp %s test_clm.nc' % clm_file_name, shell=True, stdout=PIPE)
#		sts = os.waitpid(p.pid, 0)[1]
#		clm_file_name = 'test_clm.nc'
	
	
	bry_ncfile = nc.Dataset(bry_file_name,'a')
	ini_ncfile = nc.Dataset(ini_file_name,'a')
	grid_ncfile = nc.Dataset(grid_file_name,'r')
#	clm_ncfile = nc.Dataset(clm_file_name,'a')
	
	# make variables in the boundary file
	new_var_copy_of_old_var(bry_ncfile,'salt_time','nitrate_time')
	new_var_copy_of_old_var(bry_ncfile,'salt_time','photoplankton_time')
	sides = ['north','south','east','west']
	for side in sides:
		uniform_new_var_from_old(bry_ncfile,'salt_%s' % side,'nitrate_%s' % side,fill_val=0.01,long_name='Nitrate concentration',units='umol/l')
		new_var_linear_function_of_old_var(bry_ncfile,'salt_%s' % side,'photoplankton_%s' % side, a=0.01, b=1.0/30.0, long_name='Photoplankton concentration',units='bugs/l')
	
	# make new variables in the initial conditions file
	uniform_new_var_from_old(ini_ncfile,'salt' ,'nitrate',fill_val=0.01,long_name='Nitrate concentration',units='umol/l')
	new_var_linear_function_of_old_var(ini_ncfile,'salt','photoplankton', long_name='Photoplankton concentration',units='bugs/l')
	
	
	new_var_function_of_old_var_and_depth(bry_ncfile,'salt_north','var1_north',grid_ncfile)
	new_var_function_of_old_var_and_depth(bry_ncfile,'salt_south','var1_south',grid_ncfile)
	new_var_function_of_old_var_and_depth(bry_ncfile,'salt_east','var1_east',grid_ncfile)
	new_var_function_of_old_var_and_depth(bry_ncfile,'salt_west','var1_west',grid_ncfile)
	
	new_var_function_of_old_var_and_depth(ini_ncfile,'salt','var1',grid_ncfile)
	
	# close all of the netCDF files
	bry_ncfile.close()
	ini_ncfile.close()
#	clm_ncfile.close()
	grid_ncfile.close()