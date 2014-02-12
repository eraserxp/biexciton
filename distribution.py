from numpy import *
import os

number_of_molecule = 501
time = raw_input("Enter the time (in unit of 1e-7 second): ")
time = int(time) + 1
time = str(time)
# print line number 52
# sed -n '52p' file 
os.system("sed -n '" + time+ "p'"+ " biex-250.dat | awk '{ print $2 }' > distribution.dat ")

for i in arange(-number_of_molecule/2+1, number_of_molecule/2 + 1):
	file_name2 = 'biex' + str(i) + '.dat'
	commandline = "sed -n '" + time + "p' " + file_name2 + \
 " | awk '{ print $2 }' >> distribution.dat"
	os.system(commandline) 
	print i

#normalized it
p_k = loadtxt("distribution.dat", unpack=True)
p_k = p_k/number_of_molecule
k = arange(-number_of_molecule/2, number_of_molecule/2 + 1)*pi/number_of_molecule 
savetxt('p_vs_k.dat', transpose((k, p_k)) )
os.system("xmgrace p_vs_k.dat")
