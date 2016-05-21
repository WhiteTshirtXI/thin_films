import matplotlib.pyplot as plt 
from matplotlib.pyplot import cm
import numpy as np
import glob
from scipy import interpolate, interpolate
from math import log10, floor, pi

def readFiles(path) :
	files = sorted(glob.glob(path))
	nfiles = len(files)

	if nfiles > 0 :
		for f in files :
			x,u = np.loadtxt(f, unpack=True,skiprows=2)
			plt.plot(x,u,'--.')
		plt.show()
	
#	if nfiles > 0 :
#		f = files[-1]
#		x,u = np.loadtxt(f, unpack=True,skiprows=2)
#		plt.plot(x,u)
#		plt.show()

path = "../output/*"
readFiles(path)

