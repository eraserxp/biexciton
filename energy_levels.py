import matplotlib
matplotlib.rc('text', usetex=True)
from pylab import *

ka, En2 = loadtxt('N2_dispersion.dat').T
ka, Exx = loadtxt('biexciton_spectrum.dat').T
ka, Emin = loadtxt('Emin_continuum.dat').T
ka, Emax = loadtxt('Emax_continuum.dat').T
plot(ka, En2, ka, Exx, ka, Emin, ka, Emax)
xlim(-pi,pi) 
#xlabel(r'\font\a ptmri8r at 14pt\a ' + 'ka')
xlabel(r'\font\a ptmri8r at 18pt\a ' + 'ka')
#xlabel('ka', style='italic') #Times-Italic font in xmgrace
ylabel('Energy (kHz)') 
fill_between(ka, Emin, Emax, alpha=0.5,color='yellow')
show()
