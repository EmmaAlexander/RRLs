import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import scipy.stats as stats

#datafiles = ['411to434stacked.txt','435to458stacked.txt','459to482stacked.txt','483to506stacked.txt','507to530stacked.txt']
#markers = ['ob','Dg','^r','*c','pm']
datafiles = ['411to422stacked10kms.txt','423to434stacked10kms.txt','435to446stacked10kms.txt','447to458stacked10kms.txt','459to470stacked10kms.txt','471to482stacked10kms.txt','483to494stacked10kms.txt','495to506stacked10kms.txt','507to518stacked10kms.txt','519to530stacked10kms.txt','531to550stacked10kms.txt']
#markers = ['om']

#datafiles = ['']
outfilename = 'datafits10kmsalpha.txt'

numfiles = len(datafiles)

minchans = np.empty(shape=(numfiles))
maxchans = np.empty(shape=(numfiles))
x0s = np.empty(shape=(numfiles))
amps = np.empty(shape=(numfiles))
fwhm_Gs = np.empty(shape=(numfiles))
fwhm_Ls = np.empty(shape=(numfiles))
fwhms = np.empty(shape=(numfiles))
redchisquared = np.empty(shape=(numfiles))
pltx = np.arange(start=-200,stop=200,step=0.5)

outfile = open(outfilename,'w') 
outfile.write('#Fitting parameters for line profiles\n')
outfile.write('#Min transition number; max transition number, x0, amplitude, fhwm G, fwhm L, fwhm, reduced chi squared\n')

for i in range(0,numfiles):
	filename = datafiles[i]
	mintran = filename[0:3]
	maxtran = filename[5:8]
	plotlabel = 'n = ' + mintran + ' to ' + maxtran
	data = np.loadtxt(filename)
	x = data[:,0]
	y = -data[:,1]
	yerrs = data[:,2]

	# Fit the data using a Gaussian
	g_init = models.Gaussian1D(amplitude=1, mean=13, stddev=20.)
	fit_g = fitting.LevMarLSQFitter()
	g = fit_g(g_init, x, y)
	print g

	# Fit the data using a Voigt
	#v_init = models.Voigt1D(x_0=10.,amplitude_L=0.005,fwhm_L=25.,fwhm_G=20.)
	v_init = models.Voigt1D(x_0=10.,amplitude_L=0.005,fwhm_L=1.,fwhm_G=25.)
	fit_v = fitting.LevMarLSQFitter()
	v = fit_v(v_init, x, y)

	x0s[i]= v.x_0.value
	amps[i] = v.amplitude_L.value
	fwhm_Gs[i] = v.fwhm_G.value
	fwhm_Ls[i] = v.fwhm_L.value

	fwhms[i] = 0.5346*fwhm_Ls[i] + np.sqrt((0.2166*fwhm_Ls[i]*fwhm_Ls[i])+(fwhm_Gs[i]*fwhm_Gs[i]))

	resids = np.subtract(v(x),y)

	chisquared, pval = stats.chisquare(y,f_exp=v(x))
	redchisquared[i] = chisquared/(len(y)-1)

	s = mintran+' '+maxtran+' '+str(x0s[i])+' '+str(amps[i])+' '+str(fwhm_Gs[i])+' '+str(fwhm_Ls[i])+' '+str(fwhms[i])+' '+str(redchisquared[i])+'\n'
	outfile.write(s)

	plt.plot(pltx,-v(pltx),label=plotlabel)
	#PLOT
	ax1 = plt.subplot2grid((4,1),(0,0),rowspan=3)
	ax1.errorbar(x,-y,yerr=yerrs,linestyle='none')
	ax1.plot(x, -y, 'om', label=plotlabel)
	ax1.plot(pltx, -v(pltx),color='b')
	plt.xlim(np.nanmin(x)-5,np.nanmax(x)+5)
	#plt.ylim(-0.001,0.0007)
	plt.ylabel('Optical depth')
	plt.axhline(y=0,color='k')
	plt.axvline(x=0,color='k')

	ax2 = plt.subplot2grid((4,1),(3,0))
	ax2.errorbar(x, resids,yerr=yerrs,linestyle='none')
	ax2.plot(x, resids, 'om')
	#plt.xlim(-200,200)
	plt.xlim(np.nanmin(x)-5,np.nanmax(x)+5)
	plt.ylim(-0.0006,0.0006)
	plt.xlabel(r'$V_{LSR}$ [km/s]')
	plt.ylabel('Residual')
	plt.axhline(y=0,color='k')
	plt.axvline(x=0,color='k')

	plotname = str(mintran) + 'to' + str(maxtran) + 'alpha10kms.png'
	plt.legend()
	plt.savefig(plotname,dpi=500)
	plt.show()

plt.xlabel('Velocity [km/s]')
plt.ylabel('Optical depth')
plt.legend()
plt.show()

# Two subplots, the axes array is 1-d
#f, axarr = plt.subplots(2, sharex=True)
#axarr[0].plot(x, y)
##axarr[0].set_title('Sharing X axis')
#axarr[1].scatter(x, y)

#plt.plot(redchisquared)
#plt.show()
