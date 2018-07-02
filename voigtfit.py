import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting

data1 = np.loadtxt('50to60MHz.txt')

x = data1[:,0]
y = -data1[:,1]

# Fit the data using a Gaussian
g_init = models.Gaussian1D(amplitude=1, mean=13, stddev=20.)
fit_g = fitting.LevMarLSQFitter()
g = fit_g(g_init, x, y)
print g

# Fit the data using a Gaussian
v_init = models.Voigt1D(x_0=10.,amplitude_L=0.005,fwhm_L=25.,fwhm_G=20.)
fit_v = fitting.LevMarLSQFitter()
v = fit_v(v_init, x, y)
print v

# Plot the data with the best-fit model
plt.plot(x, -y, 'ko')
plt.plot(x, -g(x), label='Gaussian')
plt.plot(x, -v(x), label='Voigt',color='r')
plt.xlim(-100,100)
plt.xlabel('Velocity [km/s]')
plt.ylabel('Optical depth')
plt.legend(loc=2)

plt.show()