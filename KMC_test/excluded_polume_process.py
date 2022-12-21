from KMCLib import *
import numpy as np

from KMCLib import *
# Setup a diffusion process up-down, middle
#  *,  E,  *  ->   *,  P1, *
#  Pi, Pi+1, Pi+2 ->   P0, E,  P2
#  *,  *,  *  ->   *,  *,  *
#                   P0             P1            P2           Moved
coordinates = [[0.0, 0.0, 0.0],[2.0, 0.0, 0.0],[4.0,0.0,0.0],[1.0,1.0,0.0]]
p0 = KMCProcess(coordinates=coordinates,
                elements_before=['P0','P1','P2','E'],
                elements_after=['P0','E','P2','P1'],
                move_vectors=None,
                basis_sites=[0],
                rate_constant=1.)
