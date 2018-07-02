import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting

#datafiles = ['411to434stacked.txt','435to458stacked.txt','459to482stacked.txt','483to506stacked.txt','507to530stacked.txt']
#markers = ['ob','Dg','^r','*c','pm']
datafiles = ['507to530stacked.txt']
markers = ['pm']


for i in range(0,len(datafiles)):
	filename = datafiles[i]
	mintran = filename[0:3]
	maxtran = filename[5:8]
	plotlabel = 'n = ' + mintran + ' to ' + maxtran
	data = np.loadtxt(filename)
	x = data[:,0]
	y = -data[:,1]

	# Fit the data using a Gaussian
	#g_init = models.Gaussian1D(amplitude=1, mean=13, stddev=20.)
	#fit_g = fitting.LevMarLSQFitter()
	#g = fit_g(g_init, x, y)
	#print g

	# Fit the data using a Voigt
	#v_init = models.Voigt1D(x_0=10.,amplitude_L=0.005,fwhm_L=25.,fwhm_G=20.)
	v_init = models.Voigt1D(x_0=10.,amplitude_L=0.005,fwhm_L=1.,fwhm_G=25.)
	fit_v = fitting.LevMarLSQFitter()
	v = fit_v(v_init, x, y)

	# Plot the data with the best-fit model
	#plt.plot(x, -y, markers[i], label=plotlabel)
	#plt.plot(x, -g(x), label='Gaussian')
	#plt.plot(x, -v(x))

	ax1 = plt.subplot2grid((3,1),(0,0),rowspan=2)
	ax1.plot(x, -y, markers[i], label=plotlabel)
	ax1.plot(x, -v(x))

	ax2 = plt.subplot2grid((3,1),(2,0))
	ax2.plot(x, -v(x))

plt.xlim(-100,120)
plt.xlabel('Velocity [km/s]')
plt.ylabel('Optical depth')
plt.axhline(y=0,color='k')
plt.axvline(x=0,color='k')
plt.legend()

plt.show()

# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(x, y)
axarr[0].set_title('Sharing X axis')
axarr[1].scatter(x, y)
