#!/usr/bin/env python
from subprocess import Popen, PIPE
from numpy import *
import re

""" find the distance at the minimum energy """

# function to compute bonding energy for a diatomic system at a given distance
def EnergyCal(z1, z2, dist):
    #uses the subprocess module to run the external program 'diatomic_energy.py'
    proc = Popen('python diatomic_energy.py', stdin = PIPE, stdout = PIPE, shell = True)
    output = proc.communicate('%i %i %f' %(z1, z2, dist))[0]    # pass three variables to stdin

    # use a regular expression to obtain the value of energy
    label = r'Potential Energy:\s*([-]?\d*\.?\d*)'
    Energy = re.search(label, output, re.MULTILINE).group(1)

    return (float(Energy))


# function to find the optimal distance
def OptimalDistance(z1, z2):
    # create an  array of distances, range from 0.5 to 10, by step 0.001
    d = arange(0.5, 10, 0.5)
    
    # initialize the minimum energy and distance
    minEnergy = EnergyCal(z1,z2,d[0])
    optDist = d[0]

    # loop over array to find the minimum energy and the corresponding optimal distance
    for i in range (1, len(d)):
        if EnergyCal(z1,z2,d[i]) < minEnergy:
            minEnergy = EnergyCal(z1,z2,d[i])
            optDist = d[i]

    d2 = arange(optDist-0.5, optDist+0.5, 0.01)
    for i in range (1, len(d2)):
        if EnergyCal(z1,z2,d2[i]) < minEnergy:
            minEnergy = EnergyCal(z1,z2,d2[i])
            optDist = d2[i]

    d3 = arange(optDist-0.01, optDist+0.01, 0.001)
    for i in range (1, len(d3)):
        if EnergyCal(z1,z2,d3[i]) < minEnergy:
            minEnergy = EnergyCal(z1,z2,d3[i])
            optDist = d3[i]
    return optDist


""" main script """
atomicNum1 = 1
atomicNum2 = 2
distance = 1.3

print 'The minimum energy: %f kcal/mol' %(EnergyCal(atomicNum1, atomicNum2, distance))
print 'The optimal distance:', OptimalDistance(atomicNum1, atomicNum2)
