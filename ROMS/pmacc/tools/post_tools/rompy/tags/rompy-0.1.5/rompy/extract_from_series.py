import os
import datetime as dt
import time
import glob

import pytz
import netCDF4 as nc
import numpy as np

import extract_utils
import extract_from_file

__version__ = '0.1.5'

# def file_time(f):
# 	UTC = pytz.timezone('UTC')
# 	ncf = nc.Dataset(f,mode='r')
# 	ot = ncf.variables['ocean_time']
# 	base_time = dt.datetime.strptime(ot.units,'seconds since %Y-%m-%d %H:%M:%S').replace(tzinfo=UTC)
# 	offset = dt.timedelta(seconds=ot[0][0])
# 	ncf.close()
# 	return base_time + offset
# 
# def file_timestamp(f):
# 	t = file_time(f)
# 	return time.mktime(t.timetuple())
# 
# def calc_num_time_slices(file1,file2,td):
# 	d1 = file_time(file1)
# 	d2 = file_time(file2)
# 	
# 	gap = d2-d1
# 	gap_seconds = gap.days*60*60*24 + gap.seconds
# 	td_seconds = td.days*60*60*24 + td.seconds
# 	
# 	t1 = time.mktime(d1.timetuple())
# 	t2 = time.mktime(d2.timetuple())
# 	tgap = t2 - t1
# 	
# 	print(gap_seconds, tgap)
# 	
# 	return int(gap_seconds/td_seconds + 1)
# 
# def filelist_from_datelist(datelist, basedir='/Users/lederer/Repositories/PSVS/rompy/', basename='ocean_his_*.nc'):
# 	files = glob.glob(os.path.join(basedir,basename))
# 	files.sort()
# 	
# 	master_file_list = []
# 	master_timestamp_list = []
# 	master_datetime_list = []
# 	timelist = []
# 	returnlist = []
# 	
# 	for file in files:
# 		master_file_list.append(file)
# 		master_timestamp_list.append(file_timestamp(file))
# 		master_datetime_list.append(file_time(file))
# 	
# 	tsarray = np.array(master_timestamp_list)
# 
# 	for d in datelist:
# 		try:
# 			t = time.mktime(d.timetuple())
# 			timelist.append(t)
# 			tl = np.nonzero(tsarray <= t)[0][-1]
# 			th = np.nonzero(t <= tsarray)[0][0]
# 			if not tl == th:
# 				fraction = (t-tsarray[th])/(tsarray[tl]-tsarray[th])
# 			else:
# 				fraction = 0.0
# 			
# 			returnlist.append({
# 							'file0': master_file_list[tl],
# 							'file1': master_file_list[th],
# 							'fraction': fraction,
# 							'timestamp':t,
# 							'datetime':d
# 							})
# 		except IndexError, e:
# 			print('Index out of bounds for %s' % (d.isoformat()))
# 			
# 	return returnlist

def extract_from_series(file_list,extraction_type='point',varname='zeta',**kwargs):
	UTC = pytz.timezone('UTC')
	pacific = pytz.timezone('US/Pacific')
	
	#print('Hello from extractFromSeries')
	
	if extraction_type == 'point':
		# assume var is a 2d var for now
		file_list.sort()
		
		x = kwargs['x']
		y = kwargs['y']
		
		data = np.zeros((len(file_list),len(x)))
		for i in range(len(x)):	
			data[0,i],junk = extract_from_file.extract_from_file(file_list[0],varname=varname,extraction_type='point',x=x[i],y=y[i])
		
		f1 = file_list[0]
		f2 = file_list[1]
		time_list = [extract_utils.file_time(file_list[0])]
		
		for i in range(1,len(file_list)):
			for j in range(len(x)):
				data[i,j],junk = extract_from_file.extract_from_file(file_list[i],varname=varname,extraction_type='point',x=x[j],y=y[j])
			time_list.append(extract_utils.file_time(file_list[i]))
		return (data,time_list)
		
#		for i in range(1,ntimes):
#			time_list.append(time_list[-1]+freq)
		
	
	if extraction_type == 'profile':
		
		file_list.sort()
		
		ocean_time = None
		
		data = None
		z = None
		
		x = kwargs['x']
		y = kwargs['y']
		
		if len(x) > 1:
			x = np.array(x[0])
		else:
			x = np.array(x)
		if len(y) > 1:
			y = np.array(y[0])
		else:
			y = np.array(y)

		for i in range(len(file_list)):
			file = file_list[i]
			
			if not os.path.exists(file):
				raise IOError('File %s could not be located on the filesystem' %file)
			ncf = nc.Dataset(file,mode='r')
			
			if varname not in ncf.variables:
				raise IOError('File %s does not have a variable named %s' % (file, varname))
			ncf.close()
			
			d,junk = extract_from_file.extract_from_file(file_list[i],varname=varname,extraction_type='profile',x=x,y=y)
			
			if data == None:
				data = np.zeros((d.shape[0],len(file_list)))
			data[:,i] = d.T
			
			if z == None:
				z = np.zeros(data.shape)
			z[:,i] = junk['zm'].T
			
#				ocean_time = np.zeros(data.shape)
			if ocean_time == None:
				ocean_time = np.zeros(data.shape)
			ot = extract_utils.file_timestamp(file_list[i])
#			ot = file_time(file)
			for j in range(data.shape[0]):
				ocean_time[j,i] = ot
		return (data, ocean_time,z)
	return

def extract_from_datetime_list(datelist,x,y,varname='salt',basedir='./',**kwargs):
	filelist = extract_utils.filelist_from_datelist(datelist, basedir=basedir)
	for f in filelist:
		(d0,ot0,z0) = extract_from_series([f['file0']],x=x,y=y,varname=varname,extraction_type='profile',**kwargs)
		(d1,ot1,z1) = extract_from_series([f['file1']],x=x,y=y,varname=varname,extraction_type='profile',**kwargs)
		
		di = np.zeros(d0.shape)
		oti = np.zeros(d0.shape)
		zi = np.zeros(d0.shape)

		f['d'] = d0 + (d1-d0)*f['fraction']
		f['ot'] = ot0 + (ot1-ot0)*f['fraction']
		f['z'] = z0 + (z1-z0)*f['fraction']
	
	d = np.zeros((filelist[0]['d'].shape[0], len(filelist)*filelist[0]['d'].shape[1]))
	ot = np.zeros((filelist[0]['d'].shape[0], len(filelist)*filelist[0]['d'].shape[1]))
	z = np.zeros((filelist[0]['d'].shape[0], len(filelist)*filelist[0]['d'].shape[1]))

	for i in range(len(filelist)):
		d[:,i] = filelist[i]['d'].T
		ot[:,i] = filelist[i]['ot'].T
		z[:,i] = filelist[i]['z'].T
	return (d,ot,z)

def extract_from_two_datetimes(x,y,dt0,dt1,varname='salt',interval=3600,basedir='./',**kwargs):
	t0 = time.mktime(dt0.timetuple())
	t1 = time.mktime(dt1.timetuple())
	interval = float(interval)

	time_stamps = interval*np.arange(t0/interval,(t1/interval)+1)

	date_times = []
	for ts in time_stamps:
		date_times.append(dt.datetime.fromtimestamp(ts))
#	for d in date_times:
#		print(d.isoformat())
	print(dt0.isoformat(),date_times[0].isoformat())
	print(dt1.isoformat(),date_times[-1].isoformat())
	return extract_from_datetime_list(date_times,x,y,varname=varname,basedir=basedir,**kwargs)
	