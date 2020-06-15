from sys import argv
import csv
import numpy as np
import os
from math import radians, cos, sin, asin, sqrt, exp
import copy
import random
import pyemd
from pyemd import emd
from pyemd import emd_with_flow

os.chdir('/Volumes/Marius_SSD/American-Flyway/Connectivity_NAbirds/Redistribution-model')

# Load species seasonal abundance distributions (estimated from eBird data)
abundance_BR = np.loadtxt('Data/STEMs/seasonalAbundance_wlswar_BR.csv', delimiter=';')
abundance_NB = np.loadtxt('Data/STEMs/seasonalAbundance_wlswar_NB.csv', delimiter=';')
abundance_BR = abundance_BR / sum(abundance_BR)
abundance_NB = abundance_NB / sum(abundance_NB)

# Load matrix of pairwise distance between every hexagons on the grid
distanceMatrix = np.loadtxt('ideal-optimal-redistribution/distanceMatrix.csv', delimiter=';')

# Compute optimal redistribution using the Earth Mover's Distance algorithm 
EMD_results = emd_with_flow(abundance_BR, abundance_NB, distanceMatrix)
flow=0
for i in range(0,len(distanceMatrix)):
    flow = flow + sum(EMD_results[1][i])
EMD_results2 = EMD_results[0] / flow

print EMD_results2

# Save simulated migratory connectivity 
np.savetxt("ORSIM-outputs/ORSIMresults_wlswar.csv", EMD_results[1], delimiter=',')
