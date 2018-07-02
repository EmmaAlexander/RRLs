# Reads in EDA fits files and finds radio recombination lines 
# Emma Alexander, June-Sept 2017, ASTRON

# import modules
import astropy.io.fits as fits
from astropy.time import Time
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt 
plt.rcParams['agg.path.chunksize'] = 20000 # soemthing for matplotlib
import glob
import scipy.signal as scisig
from scipy.stats import linregress
import math as mth
from astropy.modeling import models, fitting
import scipy.stats as stats

#######################################################
# Parameters that need hardcoding/ checking before use

# parameters
# taken from fits file, 
steptime = 2.0032 #time per scan in seconds
stepfreq = 0.00125 #channel bandwidth in MHz
reffreq = 0 # Frequency of lowest channel in MHz
c = 299792458 # constant 
velcorrection = 20.57 #km/s - specific for this observation only 

numfiles = 10 # number of fits files 
tperfile = 600 # number of time steps per fits file
numchans = 262144 # total number of frequency channels 
# ^(could probably determine all these automatically)

# create a list for filenames, the data, and time 
filenames = []
datalist = []
avdatalist = []
timelist = []

# specify the range of n of the RRLs to be stacked
# note: these refer to carbon ns, adaption would be needed otherwise 
# lines stacked include the min and max, with all inbetween. 
maxlinenum = 360#696#422
minlinenum = 334#519#411

# paramter to determine over how many freq chans baseline is fitted:
freqsplit = 256
# parameter to determine how much of spectrum either side of lines is sampled:
lineparam = 120
# velocity resolution of final bined stack: 
velocityres = 10 #km/s

#####################################################

# load in data for RRL lines and their frequencies: 
# Carbon alpha:
x,y,transnum,lines = np.loadtxt('RRL_CI_alpha_d_30_330_MHz_Oonk.txt', dtype='string,string,float,float',unpack=True)
#x,y,transnum,lines = np.loadtxt('RRL_CI_beta_d_30_330_MHz_Oonk.txt', dtype='string,string,float,float',unpack=True)
# Hydrogen alpha:
xb,yb,transnumb,linesb = np.loadtxt('RRL_CI_beta_d_30_330_MHz_Oonk.txt', dtype='string,string,float,float',unpack=True)
# Carbon beta:
xH,yH,transnumH,linesH = np.loadtxt('RRL_HI_alpha_d_30_330_MHz_Oonk.txt', dtype='string,string,float,float',unpack=True)

# select lines in the right range
lines = lines[(transnum >= minlinenum)]# & transnum<= maxlinenum)] # only want relevent lines
transnum = transnum[(transnum >= minlinenum)]# & transnum<= maxlinenum)]  # and their relevent transition numbers
lines = lines[(transnum <= maxlinenum)]
transnum = transnum[(transnum <= maxlinenum)]

numlines = len(lines)

avlinefreq = np.nanmean(lines) # average frequency of the line selection
avn = np.nanmean(transnum) # average n of the selection 

mintran = np.min(transnum)
maxtran = np.max(transnum)

# find the frequencies over which the n selection falls 
globminfreq = np.min(lines) - (freqsplit/2)*stepfreq 
globmaxfreq = np.max(lines) + (freqsplit/2)*stepfreq 
# round to nearest multiple of frequencysplit to find corresponding channel number
globminchan = int(freqsplit*np.floor(np.divide((np.min(lines)/stepfreq - (freqsplit/2)),freqsplit)))
globmaxchan = int(freqsplit*np.ceil(np.divide((np.max(lines)/stepfreq + (freqsplit/2)),freqsplit)))

# get relevnt HRRLS and beta CRRLs
linesH = linesH[(linesH>=globminfreq)]
linesH = linesH[linesH<=globmaxfreq]
linesb = linesb[(linesb>=globminfreq)]
linesb = linesb[linesb<=globmaxfreq]

# find the nearest channels to the lines 
nearestinds = np.rint(np.divide(lines,stepfreq)) - globminchan

# find the frequency range of the selection 
freqrange = globmaxchan-globminchan # MHz
numfreqchunks = int(freqrange/freqsplit)

# array of all the relevent frequencies:
allfreqs= np.multiply(np.linspace(globminchan,globmaxchan,freqrange),stepfreq) 

# find all the wanted files
print "Finding and reading data files"
for infile in glob.glob('specdump_*_ch0.fits'):

	# convert filename to a string
	filestr = str(infile)
	# append this to the filename list
	list.append(filenames,filestr)
	# find the timestamp of the file from the name
	time = filestr.replace('specdump_','')
	time = time.replace('_','')
	time = time.replace('ch0.fits','')
	list.append(timelist,time)

	#open and read fits file
	hdu = fits.open(infile)
	header = hdu[0].header
	data = hdu[0].data
	hdu.close()

	# store data in the list
	list.append(datalist,data[:,globminchan:globmaxchan])

print "All data read"
if len(timelist) == 0:
	print "Error: no files found"

# make empty arrays to concatenate data into 
alldata = np.empty([tperfile*numfiles])
alltimeMJD = np.empty([tperfile*numfiles])

#convert channels into frequency
allfreqs= np.multiply(np.linspace(0,numchans,numchans),stepfreq) + reffreq 

# convert time to MJD and combine into one continous 
timesteps = np.divide(np.multiply(np.linspace(0,tperfile-1,tperfile),steptime),3600.*24.)

##########################
# the following is not neccessary for the code in its current form, but could be useful later
# puts a time stamp on all the time steps 
#for i in range(0,numfiles):
#	timestr = str(timelist[i])
#	timestr2 = timestr[0:4]+'-'+timestr[4:6]+'-'+timestr[6:8]+' '+timestr[8:10]+':'+timestr[10:12]+':'+timestr[12:14]
#	t = Time(timestr2,format='iso', scale='utc')
#	alltimeMJD[i*tperfile:(i+1)*tperfile] = timesteps + t.mjd

## find time of first scan in MJD
#reftime = np.nanmin(alltimeMJD)
# convert all times into seconds since reference time
#alltime = np.multiply(np.subtract(alltimeMJD,reftime),24.*3600.)
###########################

# concatenate data from all files 
print "Concatenating data"
alldata = np.ma.concatenate(datalist)
masked_data = np.copy(alldata)

# temporarily mask line frequencies to fit bsaline 
print "Masking lines to fit baseline"
for i in range(0,numlines):
	minmaskfreq = lines[i] - 0.01 #might need to tweak 0.01
	maxmaskfreq = lines[i] + 0.01

	minmaskfreqH = linesH[i] - 0.01
	maxmaskfreqH = linesH[i] + 0.01

	minmaskchan = int(np.rint(np.divide(minmaskfreq,stepfreq))) - globminchan
	maxmaskchan = int(np.rint(np.divide(maxmaskfreq,stepfreq))) - globminchan

	minmaskchanH = int(np.rint(np.divide(minmaskfreqH,stepfreq))) - globminchan
	maxmaskchanH = int(np.rint(np.divide(maxmaskfreqH,stepfreq))) - globminchan

	masked_data[:,minmaskchan:maxmaskchan]=np.ma.masked
	masked_data[:,minmaskchanH:maxmaskchanH]=np.ma.masked

print "Normalising & averaging"

# manual selection of time steps that are okay to use (start, stop)
timechunks = ([(0,29),(31,92),(93,153),(155,214),(235,276),(278,336),(339,399),(401,461),(480,541),(542,600),(603,661),
	(665,724),(726,784),(786,847),(848,907),(909,968),(970,1030),(1032,1091),(1093,1153),(1155,1215),(1233,1291),
	(1294,1353),(1372,1431),(1433,1493),(1495,1554),(1557,1616),(1619,1678),(1681,1740),(1742,1802),(1804,1864),
	(1866,1926),(1928,1987),(1990,2049),(2052,2111),(2113,2173),(2175,2234),(2237,2296),(2298,2357),(2359,2481),
	(2483,2542),(2545,2604),(2607,2666),(2669,2699),(2701,2728),(2730,2790),(2792,2819),(2822,2851),(2854,2913),
	(2916,2975),(2978,3037),(3039,3099),(3101,3160),(3163,3222),(3224,3284),(3286,3346),(3348,3408),(3410,3469),
	(3472,3531),(3535,3593),(3596,3655),(3657,3717),(3719,3778),(3781,3840),(3844,3902),(3904,3963),(3966,4026),
	(4028,4087),(4090,4149),(4152,4210),(4213,4272),(4274,4333),(4336,4396),(4398,4458),(4460,4520),(4522,4581),
	(4584,4643),(4645,4705),(4708,4767),(4769,4828),(4830,4890),(4892,4952),(4934,5013),(5015,5074),(5077,5136),
	(5138,5197),(5199,5259)])

# empty, to fill 
averagelist = []
newfreqlist = []

masklines = np.all(masked_data,axis=0)

# loop over each frequency selection to fit 
for j in range(0,numfreqchunks):
	minchan = j*freqsplit
	maxchan = (j+1)*freqsplit
	x = np.arange(start=minchan,stop=maxchan,step=1)
	xfreqs = np.multiply(x,stepfreq) + globminfreq
	normalisedlist = []

	# loop over each time chunk 
	for i in range(0,len(timechunks)):
		tmin=timechunks[i][0]
		tmax=timechunks[i][1]
		masklineschunk = masklines[minchan:maxchan]
		# y data to fit: lines are masked
		yf = np.nanmean(masked_data[tmin:tmax,minchan:maxchan],axis=0)
		yf = np.ma.masked_where(masklineschunk!=True,yf)
		# true y data: 
		yt = np.nanmean(alldata[tmin:tmax,minchan:maxchan],axis=0)
		# remove mask from true y that somehow has snuck in?
		yt.mask = False
		# fit the baselines 
		p = np.ma.polyfit(x,yf,3)
		#f = p[0]*np.power(x,4) + p[1]*np.power(x,3)+ p[2]*np.power(x,2)+ p[3]*x +p[4]
		f = np.power(x,3)*p[0] + np.power(x,2)*p[1] + x*p[2] + p[3]
		# normalise
		n = np.divide(yt,f) - 1
		list.append(normalisedlist,n)

	concated = np.concatenate(normalisedlist)
	normalised = np.reshape(np.concatenate(normalisedlist),(len(timechunks),freqsplit))
	averaged = np.nanmean(normalised,axis=0)

	list.append(newfreqlist,allfreqs[minchan+globminchan:maxchan+globminchan])
	list.append(averagelist,averaged)


print "Plotting"

a = np.concatenate(newfreqlist)
b = np.ma.concatenate(averagelist)

b= np.ma.masked_where(b>=0.002,b)
b= np.ma.masked_where(b<=-0.002,b)

b= np.ma.masked_where((a>=145)&(a<=147),b)
b= np.ma.masked_where((a>=226.8)&(a<=227.3),b)
b= np.ma.masked_where((a>=163.83)&(a<=164),b)
b= np.ma.masked_where((a>=149.8)&(a<=150.6),b)

# standard deviation 
std = np.nanstd(b)
#print std

# rough std noise of a quiet section of the data
threshold = 0.0003 

##################################################################
# plot
plt.plot(a,b,color='b')
plt.xlim(np.nanmin(a),np.nanmax(a))
plt.ylim(-6.5*threshold,4*threshold)
plt.axhline(y=3*threshold,color='r', linestyle='-.')
plt.axhline(y=4*threshold,color='r', linestyle='-.')
plt.axhline(y=-3*threshold,color='r', linestyle='-.')
plt.axhline(y=-4*threshold,color='r', linestyle='-.')
plt.axhline(y=-5*threshold,color='r', linestyle='-.')
plt.axhline(y=-6*threshold,color='r', linestyle='-.')

plt.axvline(x=lines[0],color='g', linestyle='-.',label=r'Cn$\alpha$')
plt.axvline(x=linesH[0],color='m', linestyle='--',label=r'Hn$\alpha$')
plt.axvline(x=linesb[0],color='c', linestyle=':',label=r'Cn$\beta$')
for i in range(1,len(lines)):
	plt.axvline(x=lines[i],color='g', linestyle='-.')
	plt.axvline(x=linesH[i],color='m', linestyle='--')
	plt.axvline(x=linesb[i],color='c', linestyle=':')
#plt.axvline(x=54.9468224313194042,color='c', linestyle=':')
#plt.axvline(x=55.2131241331094316,color='c', linestyle=':')
#plt.axvline(x=55.4811494742588067,color='c', linestyle=':')
#plt.axvline(x=55.7509124227073372,color='c', linestyle=':')

plt.axhline(y=0,color='k')
plt.plot(a,b,color='b')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Optical depth')
#plt.legend(loc=4)
#plt.savefig('absorbtionlines.png',dpi=500)
#plt.show()

##################################################

# convert to velocities
# and more plotting 

selectionlist = []
velocitylist = []

chans = np.arange(start=-lineparam,stop=lineparam,step=1)

#print len(lines)
for i in range(0,len(lines)):
	linefreq = lines[i]
	centchan = int(nearestinds[i])
	section = b[centchan-lineparam:centchan+lineparam]
	frequencysec = a[centchan-lineparam:centchan+lineparam]
	shiftedfreqs = frequencysec - linefreq
	velocities = -np.divide(c*shiftedfreqs,linefreq)/1000 + velcorrection #in km/s
	list.append(selectionlist,section)
	list.append(velocitylist,velocities)
	plt.plot(velocities,section,linestyle='none',marker='o')
plt.ylabel('Optical depth')
plt.xlabel('Velocity [km/s]')
#plt.axhline(y=0,color='k')
##plt.axvline(x=0,color='k')
#plt.xlim(-100,100)
#plt.ylim(-0.002,0.002)
#plt.savefig('velshiftedalldots.png',dpi=500)
#plt.show()
plt.clf()

shifted = np.concatenate(selectionlist).reshape(len(lines),2*lineparam)
shifted = np.ma.masked_where(shifted>0.001,shifted)
shiftedvel = np.concatenate(velocitylist).reshape(len(lines),2*lineparam)

freqrange= np.nanmax(shiftedvel[0,:]) - np.nanmin(shiftedvel[0,:])

uppervel = velocityres*mth.floor(np.nanmax(shiftedvel[0,:])/velocityres)
lowervel = velocityres*mth.ceil(np.nanmin(shiftedvel[0,:])/velocityres)

numvelchunks = int((uppervel-lowervel)/velocityres)

averagedvels = np.empty(shape=(numvelchunks))
averageddepths = np.empty(shape=(numvelchunks))
stdvels = np.empty(shape=(numvelchunks))
errdepths = np.empty(shape=(numvelchunks))

# write out to a file
#outfilename = str(int(np.rint(globminfreq))) +'to' + str(int(np.rint(globmaxfreq)))+'MHz.txt'
outfilename = str(int(mintran))+'to'+str(int(maxtran))+'NEW.txt'
outfile = open(outfilename,'w') 
header1 = '#Output of RRllinestack.py \n#Number of lines stacked = '+str(numlines)+'\n'
outfile.write(header1)
header2 = '#Min freq = '+str(globminfreq)+' Max freq = '+str(globmaxfreq)+'\n#Velocity Resolution = '+str(velocityres)+'\n'
outfile.write(header2)

# write out to file
# stacked opticl depths as a function of velocity 
for i in range(0,numvelchunks):
	minvel = lowervel+(i*velocityres)
	maxvel = lowervel+((i+1)*velocityres )
	chunkvel = shiftedvel[(shiftedvel >= minvel) & (shiftedvel <= maxvel)] 
	chunkdepth = shifted[(shiftedvel >= minvel) & (shiftedvel <= maxvel)]

	averagedvels[i] = np.nanmean(chunkvel) 
	averageddepths[i] = np.nanmean(chunkdepth)
	stdvels[i] = np.nanstd(chunkvel)
	errdepths[i] = np.nanstd(chunkdepth)/2.
	s = str(averagedvels[i]) + ' ' + str(averageddepths[i])+ ' '+ str(errdepths[i])+'\n'
	outfile.write(s)

# get stdev noise from velocities <50 km/s
noise = np.nanstd(averageddepths[(averagedvels <= -50)])
#print std

standstr = '# File end. std noise = ' + str(std) 
outfile.write(standstr)
outfile.close

#######################################################################
#fit

# change variable names for obsolete reasons (sorry)
x = averagedvels
y = averageddepths
yerrs=errdepths

# smaller spaced x for plotting purposes
pltx = np.arange(start=-250,stop=250,step=0.5)

# fit carbon peak
#g_init = models.Gaussian1D(amplitude=0.0005, mean=150, stddev=20.)
g_init = models.Gaussian1D(amplitude=0.0001, mean=5, stddev=30.)
fit_g = fitting.LevMarLSQFitter()
g = fit_g(g_init, x, y)

# fit hydrogen peak
g_initH = models.Gaussian1D(amplitude=1, mean=150, stddev=20.)
fit_gH = fitting.LevMarLSQFitter()
gH = fit_gH(g_initH, x, y)

Gmean = g.mean.value # position of peak
Gamp = g.amplitude.value # amplitude of peak
Gstddev = g.stddev.value # width of peak 

# convert this width into frequency units
convertstd = ((Gstddev*1000*avlinefreq*1000000)/c)

# integrated optical depth 
Gint = Gamp*convertstd*np.sqrt(2*mth.pi)
errGint = np.sqrt(np.square(noise)+np.square((Gint*noise)/Gamp))

errGmean = (noise*Gstddev/(2.*Gmean))
errGstddev = errGmean/2.

# signal to noise ratio of detection 
Csig = Gamp/noise

# do the same for the hydrogen peak
GmeanH = gH.mean.value
GampH = gH.amplitude.value
GstddevH = gH.stddev.value
convertstdH = ((GstddevH*1000*avlinefreq*1000000)/c)
GintH = GampH*convertstdH*np.sqrt(2*mth.pi)
errGintH = np.sqrt(np.square(noise)+np.square((GintH*noise)/GampH))
errGmeanH = (noise*GstddevH/(2.*GmeanH))
errGstddevH = errGmeanH/2.
Hsig = GampH/noise

# output quantities to terminal 
print "n_min n_max noise intC err Cmean err Cwidth err Camp err Csig "
print(minlinenum,maxlinenum,noise,Gint,errGint,Gmean,errGmean,Gstddev,errGstddev,Gamp,noise,Csig)

print "intH err Hmean err Hwidth err Hamp err Hsig "
print(GintH,errGintH,GmeanH,errGmeanH,GstddevH,errGstddevH,GampH,noise,Hsig)

# residuals from fit 
Gresids = np.subtract(y,g(x))

# chi squared, doesn't work 
#Gchisquared, Gpval = stats.chisquare(y,f_exp=g(x))
#Gredchisquared = Gchisquared/(len(y)-1)


# only various plotting code below 


plt.plot(x, y,color='k',linestyle='-.',marker='o',label='Data')
plt.plot(x, y,color='k',linestyle='none',marker='o')
plt.plot(pltx,g(pltx),color='b',label='Carbon fit')
plt.plot(pltx,gH(pltx),color='r',label='Hydrogen fit')
plt.xlim(np.nanmin(x)-5,np.nanmax(x)+5)
#plt.xlim(-250,300)
plt.ylabel('Optical depth')
plt.axhline(y=0,color='k')
plt.axvline(x=0,color='k')
plt.xlabel(r'$V_{LSR}$ [km/s]')
#plt.legend(loc=3)
plt.legend(loc=2)

#ax1 = plt.subplot2grid((4,1),(0,0),rowspan=3)
#ax1.errorbar(x,-y,yerr=yerrs,linestyle='none')
#ax1.plot(x, -y,color='k')
#ax1.plot(pltx,-g(pltx),color='b',linestyle='-.',label='Carbon')
#ax1.plot(pltx,gH(pltx),color='r',linestyle='-.',label='Hydrogen')
#ax1.plot(x,g(x),color='b',label='Gauss fit')
#plt.xlim(np.nanmin(x)-5,np.nanmax(x)+5)
#plt.ylim(-0.0004,0.0008)
#plt.xlim(-150,250)
#plt.ylabel('Optical depth')
#plt.axhline(y=0,color='k')#
#lt.axvline(x=0,color='k')

#ax2 = plt.subplot2grid((4,1),(3,0))
#ax2.errorbar(x, Gresids,yerr=yerrs,linestyle='none')
#ax2.plot(x, Gresids,color='k')
#plt.xlim(np.nanmin(x)-5,np.nanmax(x)+5)
#plt.ylim(-0.0004,0.0006)
#plt.xlabel(r'$V_{LSR}$ [km/s]')
#plt.xlim(-150,250)
#plt.ylabel('Residual')
#plt.axhline(y=0,color='k')
#plt.axvline(x=0,color='k')

plotname = str(mintran) + 'to' + str(maxtran) + 'NEW.png'
plt.tight_layout(pad=2)
#plt.legend()
plt.savefig(plotname,dpi=500)
plt.show()