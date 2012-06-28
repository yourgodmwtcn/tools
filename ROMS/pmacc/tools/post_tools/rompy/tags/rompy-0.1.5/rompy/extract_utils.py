import os
import datetime as dt
import time
import glob

import pytz
import numpy as np
import netCDF4 as nc

__version__ = '0.1.5'

def file_time(f):
	UTC = pytz.timezone('UTC')
	ncf = nc.Dataset(f,mode='r')
	ot = ncf.variables['ocean_time']
	base_time = dt.datetime.strptime(ot.units,'seconds since %Y-%m-%d %H:%M:%S').replace(tzinfo=UTC)
	offset = dt.timedelta(seconds=ot[0][0])
	ncf.close()
	return base_time + offset

def file_timestamp(f):
	t = file_time(f)
	return time.mktime(t.timetuple())

def calc_num_time_slices(file1,file2,td):
	d1 = file_time(file1)
	d2 = file_time(file2)
	
	gap = d2-d1
	gap_seconds = gap.days*60*60*24 + gap.seconds
	td_seconds = td.days*60*60*24 + td.seconds
	
	t1 = time.mktime(d1.timetuple())
	t2 = time.mktime(d2.timetuple())
	tgap = t2 - t1
	
	print(gap_seconds, tgap)
	
	return int(gap_seconds/td_seconds + 1)

def filelist_from_datelist(datelist, basedir='./', basename='ocean_his_*.nc'):
	files = glob.glob(os.path.join(basedir,basename))
	files.sort()
	
	master_file_list = []
	master_timestamp_list = []
	master_datetime_list = []
	timelist = []
	returnlist = []
	
	for file in files:
		master_file_list.append(file)
		master_timestamp_list.append(file_timestamp(file))
		master_datetime_list.append(file_time(file))
	
	tsarray = np.array(master_timestamp_list)

	for d in datelist:
		try:
			t = time.mktime(d.timetuple())
			timelist.append(t)
			tl = np.nonzero(tsarray <= t)[0][-1]
			th = np.nonzero(t <= tsarray)[0][0]
			if not tl == th:
				fraction = (t-tsarray[th])/(tsarray[tl]-tsarray[th])
			else:
				fraction = 0.0
			
			returnlist.append({
							'file0': master_file_list[tl],
							'file1': master_file_list[th],
							'fraction': fraction,
							'timestamp':t,
							'datetime':d
							})
		except IndexError, e:
			print('Index out of bounds for %s' % (d.isoformat()))
			
	return returnlist
	
def run_id(f):
	ncf = nc.Dataset(f,'r')
	rid = ncf.script_file[0:-3]
	ncf.close()
	return rid
	
def run_title(f):
	ncf = nc.Dataset(f,'r')
	title = ncf.title
	ncf.close()
	return title
	