from numpy import *

number_of_molecule = 501
t, p_n2 = loadtxt('n2ex-250.dat', usecols = (0,1), unpack=True) 
p_xx = loadtxt('biex-250.dat', usecols = [1], unpack=True)

for i in arange(-number_of_molecule/2+1, number_of_molecule/2 + 1):
	file_name1 = 'n2ex' + str(i)+ '.dat'
	file_name2 = 'biex' + str(i) + '.dat'
	p_n2_tmp = loadtxt(file_name1, usecols =[1] , unpack=True) 
	p_xx_tmp = loadtxt(file_name2, usecols =[1] , unpack=True)
	p_n2 = p_n2 + p_n2_tmp
	p_xx = p_xx + p_xx_tmp 
	print i

#convert time to micron second
t = t*10**6

#normalized it
p_n2 = p_n2/number_of_molecule
p_xx = p_xx/number_of_molecule
loss = 1 - p_n2 - p_xx
savetxt('from_n2excitation_to_biex_and_loss.txt', transpose((t, p_n2, p_xx, loss)) )



t2, p2_n2 = loadtxt('n2ex0.dat', usecols = (0,1), unpack=True) 
p2_xx = loadtxt('biex0.dat', usecols = [1], unpack=True)
loss2 = 1 - p2_n2 - p2_xx


#convert time to micron second
t2 = t2*10**6

savetxt('from_n2exciton_to_biex_and_loss.txt', transpose((t2, p2_n2, p2_xx, loss2)) )
