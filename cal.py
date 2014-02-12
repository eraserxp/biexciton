from numpy import *
Ka = 5*pi/6
J = 34.4*cos(Ka/2.0)
D = 39.8
M = 12.8*sqrt(4*D**2-J**2)*cos(Ka/2.0)/D
print "M=", M
print "tau=", 1000.0/M
