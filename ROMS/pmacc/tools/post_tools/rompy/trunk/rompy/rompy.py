import extract_from_file
import extract_from_series

__version__ = '0.1.6'

def extract(files,**kwargs):
	if files.__class__ == str:
		if kwargs.has_key('x') and kwargs.has_key('y'):
			(data, coords) = extract_from_file.extract_from_file(file=files, **kwargs)
		else:
			(data, coords) = extract_from_file.extract_from_file(files, **kwargs)
	elif files.__class__ == list:
		(data, coords) = extract_from_series.extract_from_series(files, **kwargs)
	else:
		raise TypeError('wrong type sent to the extractor')
		return False
	return (data,coords)

def test():
	print('rompy.test() works')
	extract('a.nc')
	#extract(['a.nc','b.nc'])
	try:
		extract({'first':'a.nc'})
	except TypeError, e:
		print(e)
	return

if __name__ == '__main__':
	test()
