from numpy import *

number_of_molecule = 501
t, p_n2 = loadtxt('n2ex0.dat', usecols = (0,1), unpack=True) 
p_xx = loadtxt('biex0.dat', usecols = [1], unpack=True)
loss = 1 - p_n2 - p_xx


#normalized it

savetxt('from_n2exciton_to_biex_and_loss.txt', transpose((t, p_n2, p_xx, loss)) )
