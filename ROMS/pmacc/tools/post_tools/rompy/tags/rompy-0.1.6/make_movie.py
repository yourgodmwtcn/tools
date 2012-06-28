#!/usr/bin/python

# author: Martin Michel
# created: 16.05.2008

# This script is called from an AppleScript and converts a list of
# image files into an image sequence (QuickTime movie) using the QTKit.

import os
import sys
import glob
from optparse import OptionParser

from AppKit import NSImage
from QTKit import QTMovie, QTMakeTime, QTAddImageCodecType

def movie_maker(movdir,imgdir,timeval,timescale,ih):
	print('making movie %s' % ih)
	imgpaths = glob.glob(os.path.join(imgdir,'ocean_his_*_%s.png'%ih))
	imgpaths.sort()
	
	movpath = os.path.join(movdir,'%s.mov' %ih )
	
	'''I am creating a new image sequence (QuickTime movie).'''
	mov, err = QTMovie.alloc().initToWritableFile_error_(movpath, None)
	
	if mov == None:
		errmsg = 'Could not create movie file: %s' % (movpath)
		raise IOError, errmsg
		
	duration = QTMakeTime(timeval, timescale)
	# you can also use "tiff"
#		attrs = {QTAddImageCodecType: 'avc1'}#, QTAddImageCodecQuality: codecHighQuality}
	attrs = {QTAddImageCodecType: 'avc1'}
	
	for imgpath in imgpaths:
#		img = NSImage.alloc().initWithContentsOfFile_(imgpath)
		img = NSImage.alloc().initByReferencingFile_(imgpath)
		mov.addImage_forDuration_withAttributes_(img, duration, attrs)
		del(img)
	print('Writing movie %s' %ih)
	mov.updateMovieFile()
#	del(mov, err, attrs,imgpaths, duration)

def crtimagesequence(movdir, imgdir, timeval, timescale):
	img_handles = ['hood_salt','hood_temp','hood_U','main_salt','main_temp','main_U','surface_salt','surface_temp']
#	img_handles = ['main_salt','main_temp','main_U','surface_salt','surface_temp']
#	img_handles = ['hood_salt']
	for ih in img_handles:
		movie_maker(movdir,imgdir,timeval,timescale,ih)
# 		print('making movie %s' % ih)
# 		imgpaths = glob.glob(os.path.join(imgdir,'ocean_his_*_%s.png'%ih))
# 		imgpaths.sort()
# 		
# 		movpath = os.path.join(movdir,'%s.mov' %ih )
# 		
# 		'''I am creating a new image sequence (QuickTime movie).'''
# 		mov, err = QTMovie.alloc().initToWritableFile_error_(movpath, None)
# 		
# 		if mov == None:
# 			errmsg = 'Could not create movie file: %s' % (movpath)
# 			raise IOError, errmsg
# 			
# 		duration = QTMakeTime(timeval, timescale)
# 		# you can also use "tiff"
# #		attrs = {QTAddImageCodecType: 'avc1'}#, QTAddImageCodecQuality: codecHighQuality}
# 		attrs = {QTAddImageCodecType: 'avc1'}
# 		
# 		for imgpath in imgpaths:
# 			img = NSImage.alloc().initWithContentsOfFile_(imgpath)
# 			mov.addImage_forDuration_withAttributes_(img, duration, attrs)
# 		mov.updateMovieFile()

if __name__ == '__main__':
	# file path to the new QuickTime movie
	movdir = sys.argv[1].decode('utf-8')
	
	imgdir = sys.argv[2].decode('utf-8')
	
	# time value & scale to be used for creating a QTTime duration
	timeval = float(sys.argv[3])
	
	timescale = float(sys.argv[4])
	
	crtimagesequence(movdir, imgdir, timeval, timescale)
	