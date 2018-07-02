import numpy as np
import matplotlib.pyplot as plt 

data = np.loadtxt('normalised_spectrum.txt')
freqs = data[:,0]
spec = data[:,1]

x,y,transnum,linefreqs = np.loadtxt('RRL_CI_alpha_d_30_330_MHz_Oonk.txt', dtype='string,string,float,float',unpack=True)
x,y,Htransnum,Hlinefreqs = np.loadtxt('RRL_HI_alpha_d_30_330_MHz_Oonk.txt', dtype='string,string,float,float',unpack=True)

for i in range(411,551):
	linefreq = linefreqs[(transnum == i)]
	Hlinefreq = Hlinefreqs[(transnum == i)]
	specselect = spec[(freqs<=linefreq+0.1)]
	freqselect = freqs[(freqs<=linefreq+0.1)]
	specselect = specselect[(freqselect>=linefreq-0.1)]
	freqselect = freqselect[(freqselect>=linefreq-0.1)]

	txtstring = str(int(i))+'_alpha_line.txt'
	pngstring = str(int(i))+'_alpha_line.png'

	output = np.column_stack((freqselect,specselect))
	np.savetxt(txtstring,output,fmt='%1.5f %1.8f' )

	plt.plot(freqselect,specselect)
	plt.axhline(y=0,color='k')
	plt.axvline(x=linefreq,color='g', linestyle='-.',label='C')
	plt.axvline(x=Hlinefreq,color='r', linestyle=':',label='H')
	plt.xlabel('Frequency [MHz]')
	plt.ylabel('Optical depth')
	plt.legend(loc=2)
	plt.gcf().subplots_adjust(left=0.15)
	plt.savefig(pngstring,dpi=400)
	plt.clf()
	#plt.show()