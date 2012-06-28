#!/usr/bin/env python
import os
import glob
from optparse import OptionParser

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure


from rompy import rompy, plot_utils, utils, extract_utils

def surface_map(file,img_file=None,varname='salt',clim=None):
	(data, coords) = rompy.extract(file,varname=varname,extraction_type='surface')
#	plot_utils.plot_surface(coords['xm'],coords['ym'],data)
	
	
	title = '%s %s %s %s' % ( extract_utils.run_title(file), os.path.basename(file), var_title_map[var], extract_utils.file_time(file).strftime(title_time_fmt) )
	plot_utils.plot_map(coords['xm'],coords['ym'],data,filename=img_file, clim=clim, title=title, caxis_label=clabel_map[varname])

def main_basin_curtain(file,img_file,varname,n=4,clim=None): # Main Basin
	if varname == 'U':
		main_basin_U_curtain(file,img_file,n,clim)
	else:
		x,y = utils.high_res_main_basin_xy(n=n)
		
		(data, coords) = rompy.extract(file, varname=varname, 	extraction_type='profile', x=x, y=y)
		
		title = '%s %s Main Basin %s %s' % (extract_utils.run_title(file), os.path.basename(file), var_title_map[var], extract_utils.file_time(file).strftime(title_time_fmt))
		
		plot_utils.plot_parker(coords=coords, data=data, varname=varname, 	region='Main Basin', filename=img_file, n=n, x_axis_offset=utils.offset_region(coords), clim=clim,cmap='banas_hsv_cm',labeled_contour_gap=2, title=title, caxis_label=clabel_map[varname])
	

def hood_canal_curtain(file,img_file,varname,n=1,clim=None): # Hood Canal
	if varname == 'U':
		hood_canal_U_curtain(file,img_file,n,clim)
	else:
		x,y = utils.high_res_hood_canal_xy(n=n)
		
		(data, coords) = rompy.extract(file, varname=varname, extraction_type='profile', x=x, y=y)
		
		title = '%s %s Hood Canal %s %s' % (extract_utils.run_title(file), os.path.basename(file), var_title_map[var], extract_utils.file_time(file).strftime(title_time_fmt))
		
		plot_utils.plot_parker(coords=coords, data=data, varname=varname, region='Hood Canal', filename=img_file, n=n,  x_axis_offset=utils.offset_region(coords), clim=clim, cmap='banas_hsv_cm',labeled_contour_gap=2, title=title, caxis_label=clabel_map[varname])

def hood_canal_U_curtain(file,img_file,n=1,clim=None): # velocity in Hood Canal
	x,y = utils.high_res_hood_canal_xy(n=n)
	(u, coords) = rompy.extract(file,varname='u',extraction_type='profile',x=x,y=y)
	(v, coords) = rompy.extract(file,varname='v',extraction_type='profile',x=x,y=y)
	data = np.zeros(u.shape)

	for i in range(u.shape[1]):
		if i == u.shape[1]-1:
			x_vec = np.array([x[i] - x[i-1], y[i] - y[i-1]])
		else:
			x_vec = np.array([x[i+1] - x[i], y[i+1] - y[i]])
		for j in range(u.shape[0]):
			u_vec = np.array([u[j,i], v[j,i]])
			data[j,i] = np.dot(x_vec,u_vec)/(np.sqrt(np.dot(x_vec,x_vec)))
	
	data = np.ma.array(data, mask=np.abs(data) > 100)
	
	title = '%s %s Hood Canal %s %s' % (extract_utils.run_title(file), os.path.basename(file), var_title_map['U'], extract_utils.file_time(file).strftime(title_time_fmt))
	
	hood_U_clim = (np.array(clim)/2.0).tolist()
	
	plot_utils.plot_parker(coords=coords,data=data,varname='U', region='Hood Canal', filename=img_file, n=n, clim=clim, x_axis_offset=utils.offset_region(coords), cmap='red_blue', title=title, caxis_label=clabel_map['U'])

def main_basin_U_curtain(file,img_file,n=1,clim=None): # velocity in Main Basin
	x,y = utils.high_res_main_basin_xy(n=n)
	(u, coords) = rompy.extract(file,varname='u',extraction_type='profile',x=x,y=y)
	(v, coords) = rompy.extract(file,varname='v',extraction_type='profile',x=x,y=y)
	data = np.zeros(u.shape)

	for i in range(u.shape[1]):
		if i == u.shape[1]-1:
			x_vec = np.array([x[i] - x[i-1], y[i] - y[i-1]])
		else:
			x_vec = np.array([x[i+1] - x[i], y[i+1] - y[i]])
		for j in range(u.shape[0]):
			u_vec = np.array([u[j,i], v[j,i]])
			data[j,i] = np.dot(x_vec,u_vec)/(np.sqrt(np.dot(x_vec,x_vec)))
	
	data = np.ma.array(data, mask=np.abs(data) > 100)
	
	title = '%s %s Main Basin %s %s' % (extract_utils.run_title(file), os.path.basename(file), var_title_map['U'], extract_utils.file_time(file).strftime(title_time_fmt))
	
	plot_utils.plot_parker(coords=coords,data=data,varname='U', region=' Main Basin', filename=img_file, n=n, clim=clim, x_axis_offset=utils.offset_region(coords),cmap='red_blue', title=title, caxis_label=clabel_map['U'])

def daves_curtain(file,img_file,section,varname,clim=None):
	if varname == 'U':
		daves_U_curtain(file,img_file,section,varname,clim)
	else:
		x = utils.get_daves_section_var(section=section,var='lon')
		y = utils.get_daves_section_var(section=section,var='lat')
		
		(data, coords) = rompy.extract(file, varname=varname, extraction_type='profile', x=x, y=y)
		
		title = '%s %s %s %s %s' % (extract_utils.run_title(file), os.path.basename(file), section, var_title_map[var], extract_utils.file_time(file).strftime(title_time_fmt))
		
		plot_utils.plot_parker(
			coords=coords,
			data=data,
			filename=img_file,
			title=title,
			x_axis_offset=utils.offset_region(coords),
			clim=clim,
			cmap='banas_hsv_cm',
			labeled_contour_gap=2,
			caxis_label=clabel_map[varname],
			inset=inset_dict[section],
			ctd_ind=ctd_ind_dict[section],
			label=utils.get_daves_section_var(section=section,var='label'),
			label_ind=utils.get_daves_section_var(section=section,var='label_ind')
			)
	return

def daves_U_curtain(file,img_file,section,varname,clim):
	x = utils.get_daves_section_var(section=section,var='lon')
	y = utils.get_daves_section_var(section=section,var='lat')
	(u, coords) = rompy.extract(file,varname='u',extraction_type='profile',x=x,y=y)
	(v, coords) = rompy.extract(file,varname='v',extraction_type='profile',x=x,y=y)
	data = np.zeros(u.shape)
	for i in range(u.shape[1]):
		if i == u.shape[1]-1:
			x_vec = np.array([x[i] - x[i-1], y[i] - y[i-1]])
		else:
			x_vec = np.array([x[i+1] - x[i], y[i+1] - y[i]])
		for j in range(u.shape[0]):
			u_vec = np.array([u[j,i], v[j,i]])
			data[j,i] = np.dot(x_vec,u_vec)/(np.sqrt(np.dot(x_vec,x_vec)))
	
	data = np.ma.array(data, mask=np.abs(data) > 100)
	
	title = '%s %s %s %s %s' % (extract_utils.run_title(file), os.path.basename(file), section, var_title_map['U'], extract_utils.file_time(file).strftime(title_time_fmt))

	plot_utils.plot_parker(
		coords=coords,
		data=data,
		filename=img_file,
		title=title,
		x_axis_offset=utils.offset_region(coords),
		clim=clim,
		cmap='red_blue',
		caxis_label=clabel_map[varname],
		inset=inset_dict[section],
		ctd_ind=ctd_ind_dict[section],
		label=utils.get_daves_section_var(section=section,var='label'),
		label_ind=utils.get_daves_section_var(section=section,var='label_ind')
		)

	return

# begin actual code that runs.
parser = OptionParser()

parser.add_option('-i', '--img_dir',
					dest='img_dir',
					default='./image_sequence',
					help='Location to save images. Default is ./image_sequnce')

(options, args) = parser.parse_args()

if args == []:
	fl = glob.glob('ocean_his*.nc')
	print(fl)
	file_list = [fl[0]]
else:
	file_list = args

img_dir = options.img_dir

var_list = ['salt','temp','U']

var_title_map = {'salt':'Salinity','temp':'Temperature','U':'Velocity'}
title_time_fmt = '%Y-%m-%d %H:%M UTC'

clims = {'salt':[0, 21,33, 33], 'temp': [8, 20], 'U':[-2,2]}

clabel_map = {'temp': u'\u00B0 C', 'salt': 'psu', 'U': 'm/s'}

inset_dict = {'AI_HC':'Puget Sound','AI_WB':'Puget Sound','JdF_SoG':'Strait of Georgia','JdF_SS':'Puget Sound','JdF_PS':'JdF_PS'}

ctd_ind_dict = {
		'AI_HC':utils.get_daves_section_var(section='AI_HC',var='PRISM_ctd_ind'),
		'AI_WB':utils.get_daves_section_var(section='AI_WB',var='PRISM_ctd_ind'),
		'JdF_SoG':utils.get_daves_section_var(section='JdF_SoG',var='IOS_ctd_ind'),
		'JdF_SS':utils.get_daves_section_var(section='JdF_SS',var='PRISM_ctd_ind'),
		'JdF_PS':utils.get_daves_section_var(section='JdF_PS',var='PRISM_ctd_ind')
		}
ctd_ind_dict['JdF_SoG'].extend(utils.get_daves_section_var(section='JdF_SoG',var='JEMS_ctd_ind'))

for file in file_list:
	ncf_index = os.path.basename(file)[:-3]
	print('%s' %ncf_index)
	if not os.path.exists(img_dir):	
		os.makedirs(img_dir)
	
	for var in var_list:
		hood_img_file = '%s/%s_hood_%s.png' %(img_dir, ncf_index,var)
		main_img_file = '%s/%s_main_%s.png' %(img_dir, ncf_index,var)
		surface_img_file = '%s/%s_surface_%s.png' % (img_dir, ncf_index, var)
		AI_HC_img_file = '%s/AI_HC_%s_%s.png' %(img_dir,var,ncf_index)
		AI_WB_img_file = '%s/AI_WB_%s_%s.png' %(img_dir,var,ncf_index)
		JdF_SoG_img_file = '%s/JdF_SoG_%s_%s.png' %(img_dir,var,ncf_index)
		JdF_SS_img_file = '%s/JdF_SS_%s_%s.png' %(img_dir,var,ncf_index)
		JdF_PS_img_file = '%s/JdF_PS_%s_%s.png' %(img_dir,var,ncf_index)
		
		print('making hood canal %s' % var)
		hood_canal_curtain(file, hood_img_file, var, n=8, clim=clims[var])
		print('making main basin %s' % var)
		main_basin_curtain(file, main_img_file, var, n=8, clim=clims[var])
		
		print('making AI_HC %s' % var)
		daves_curtain(file,AI_HC_img_file,section='AI_HC',varname=var,clim=clims[var])
		
		print('making AI_WB %s' % var)
		daves_curtain(file,AI_WB_img_file,section='AI_WB',varname=var,clim=clims[var])
		
		print('making JdF_SoG %s' % var)
		daves_curtain(file,JdF_SoG_img_file,section='JdF_SoG',varname=var,clim=clims[var])
		
		print('making JdF_SS %s' % var)
		daves_curtain(file,JdF_SS_img_file,section='JdF_SS',varname=var,clim=clims[var])
		
		print('making JdF_PS %s' % var)
		daves_curtain(file,JdF_PS_img_file,section='JdF_PS',varname=var,clim=clims[var])

		if not var == 'U':
			print('making surface %s' % var)
			surface_map(file,surface_img_file,var,clim=clims[var])
