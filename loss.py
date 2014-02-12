from numpy import *

number_of_molecule = 501
t, p_n2, p_xx = loadtxt('from_n2_to_biex.dat', unpack=True) 
#p_xx = loadtxt('biex0.dat', usecols = [1], unpack=True)
loss = 1 - p_n2 - p_xx


#normalized it

savetxt('loss.dat', transpose((t,  loss)) )
