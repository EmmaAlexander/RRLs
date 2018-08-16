# Reads in EDA fits files and finds radio recombination lines 
# Emma Alexander, June-Sept 2017, ASTRON
# Updated August 2018

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
import argparse, ConfigParser
c = 299792458 # constant 
maskval = 0.005

def main(args,cfg):
	# create a list for filenames, the data, and time 
	filenames = []
	datalist = []
	avdatalist = []
	timelist = []
	# empty numpy arrays
	linelist = np.empty(shape=(99999999),dtype='float')
	specieslist = np.empty(shape=(99999999),dtype=np.dtype('a16'))

	#get arguments from config file 
	steptime = float(cfg.get('input','steptime'))
	stepfreq = float(cfg.get('input','stepfreq'))
	reffreq = float(cfg.get('input','reffreq'))
	velcorrection = float(cfg.get('input','velcorrection'))
	numfiles = int(cfg.get('input','numfiles'))
	tperfile = int(cfg.get('input','tperfile'))
	numchans = int(cfg.get('input','numchans'))
	maxfreq=float(cfg.get('input','maxfreq'))
	minfreq=float(cfg.get('input','minfreq'))
	freqsplit=int(cfg.get('input','freqsplit'))
	lineparam = int(cfg.get('input','lineparam'))
	velocityres = float(cfg.get('input','velocityres'))
	speciestostack = cfg.get('input','speciestostack')

	# calculate things based on inputs
	globminchan = int(freqsplit*np.floor(np.divide((minfreq/stepfreq - (freqsplit/2)),freqsplit)))
	globmaxchan = int(freqsplit*np.ceil(np.divide((maxfreq/stepfreq + (freqsplit/2)),freqsplit)))

	# Load in the line frequencies
	totlines=0
	for line in np.genfromtxt(args.line_file_list,dtype='str'):
		species = line[0]
		filename = line[1]
		x,y,n,f = np.loadtxt(filename, dtype='string,string,float,float',unpack=True)
		for linefreq in f:
			specieslist[totlines]=str(species)
			#print specieslist[totlines]
			linelist[totlines]=linefreq
			totlines+=1 # add one to count
	# Remove empty parts of numpy arrays
	linelist=linelist[0:totlines]
	specieslist=specieslist[0:totlines]
	#print totlines

	# remove lines outside frequency range of interest 
	specieslist = specieslist[(linelist>=minfreq)]
	linelist = linelist[(linelist>=minfreq)]
	specieslist = specieslist[(linelist<=maxfreq)]
	linelist = linelist[(linelist<=maxfreq)]

	totlinesinfreq=len(linelist)

	nearestinds = np.rint(np.divide(linelist,stepfreq)) - globminchan


	# find the frequency range of the selection 
	freqrange = globmaxchan-globminchan # MHz
	numfreqchunks = int(freqrange/freqsplit)

	# Load in the .fits files
	print('Finding and reading data files')
	path = cfg.get('input','datapath')
	path= path[1:-1] #remove ' '
	for infile in glob.glob(path):

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

	print('All data read')
	if len(timelist) == 0:
		print('Error: no files found')

	# make empty arrays to concatenate data into 
	alldata = np.empty([tperfile*numfiles])
	alltimeMJD = np.empty([tperfile*numfiles])

	#convert channels into frequency
	allfreqs= np.multiply(np.linspace(0,numchans,numchans),stepfreq) + reffreq
	freqlist=allfreqs[globminchan:globmaxchan] 

	# convert time to MJD and combine into one continous 
	timesteps = np.divide(np.multiply(np.linspace(0,tperfile-1,tperfile),steptime),3600.*24.)

	print('Concatenating data')
	alldata = np.ma.concatenate(datalist)
	masked_data = np.copy(alldata)

	# find the selection of frequency and apply it to the data 

	averaged=np.nanmean(alldata,axis=0)
	#plt.plot(freqlist,averaged)

	print('Masking lines to fit baseline')
	for i in range(0,totlinesinfreq):
		minmaskfreq = linelist[i] - maskval
		maxmaskfreq = linelist[i] + maskval
		minmaskchan = int(np.rint(np.divide(minmaskfreq,stepfreq))) - globminchan
		maxmaskchan = int(np.rint(np.divide(maxmaskfreq,stepfreq))) - globminchan
		masked_data[:,minmaskchan:maxmaskchan]=np.ma.masked
		#print specieslist[i]

		#if specieslist[i] == 'Calpha':
		#	plt.axvline(x=linelist[i],color='k',label=r'Cn$\alpha$')
		#elif specieslist[i] == 'Cbeta':
		#	plt.axvline(x=linelist[i],color='b',label=r'Cn$\beta$')
		#elif specieslist[i] == 'Halpha':
		#	plt.axvline(x=linelist[i],color='g',label=r'Hn$\alpha$')
		#else:
		#	plt.axvline(x=linelist[i],color='r',label=r'other')

	#plt.xlabel('Frequency [MHz]')
	#plt.ylabel('Optical depth')
	#plt.legend(loc=4)
	#plt.show()

	print('Normalising & averaging')

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
		xfreqs = np.multiply(x,stepfreq) + minfreq
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

	newfreqs = np.concatenate(newfreqlist)
	newavs = np.ma.concatenate(averagelist)

	newavs= np.ma.masked_where(newavs>=0.002,newavs)
	newavs= np.ma.masked_where(newavs<=-0.002,newavs)

	newavs= np.ma.masked_where((newfreqs>=145)&(newfreqs<=147),newavs)
	newavs= np.ma.masked_where((newfreqs>=226.8)&(newfreqs<=227.3),newavs)
	newavs= np.ma.masked_where((newfreqs>=163.83)&(newfreqs<=164),newavs)
	newavs= np.ma.masked_where((newfreqs>=149.8)&(newfreqs<=150.6),newavs)

	# standard deviation 
	std = np.nanstd(newavs)
	#print std

	# rough std noise of a quiet section of the data
	threshold = 0.0003 

	print('Plotting')

	#plt.plot(newfreqs,newavs,color='b')
	#plt.show()

	# STACK THE LINES
	selectionlist = []
	velocitylist = []

	chans = np.arange(start=-lineparam,stop=lineparam,step=1)

	linecount = 0
	for i in range(0,totlinesinfreq):
		if specieslist[i] == speciestostack[1:-1]:
			linecount+=1
			linefreq = linelist[i]
			centchan = int(nearestinds[i])
			section = newavs[centchan-lineparam:centchan+lineparam]
			frequencysec = newfreqs[centchan-lineparam:centchan+lineparam]
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

	shifted = np.concatenate(selectionlist).reshape(linecount,2*lineparam)
	shifted = np.ma.masked_where(shifted>0.001,shifted)
	shiftedvel = np.concatenate(velocitylist).reshape(linecount,2*lineparam)

	freqrange= np.nanmax(shiftedvel[0,:]) - np.nanmin(shiftedvel[0,:])

	uppervel = velocityres*mth.floor(np.nanmax(shiftedvel[0,:])/velocityres)
	lowervel = velocityres*mth.ceil(np.nanmin(shiftedvel[0,:])/velocityres)

	numvelchunks = int((uppervel-lowervel)/velocityres)

	averagedvels = np.empty(shape=(numvelchunks))
	averageddepths = np.empty(shape=(numvelchunks))
	stdvels = np.empty(shape=(numvelchunks))
	errdepths = np.empty(shape=(numvelchunks))

	# write out to a files
	#outfilename = str(int(np.rint(globminfreq))) +'to' + str(int(np.rint(globmaxfreq)))+'MHz.txt'
	outfilename = str(int(minfreq))+'to'+str(int(maxfreq))+'MHz'+speciestostack[1:-1]+'.txt'
	outfile = open(outfilename,'w') 
	header1 = '#Output of RRllinestack.py \n#Number of lines stacked = '+str(linecount)+'\n'
	outfile.write(header1)
	header2 = '#Min freq = '+str(minfreq)+' Max freq = '+str(maxfreq)+'\n#Velocity Resolution = '+str(velocityres)+'\n'
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

	plt.plot(averagedvels,averageddepths)
	plt.show()


ap = argparse.ArgumentParser()
ap.add_argument('config_file',help='Input configuration file')
ap.add_argument('-l','--line_file_list',help='Name of text file containing line files [default linefiles.txt]',default='linefiles.txt')
args = ap.parse_args()

cfg = ConfigParser.RawConfigParser()
cfg.read(args.config_file)

main(args,cfg)
