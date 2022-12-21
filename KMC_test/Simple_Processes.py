from KMCLib import *

# Setup a diffusion process to the left.
#                   center site      second site on -1
coordinates_p0 = [[0.0, 0.0, 0.0],[0.0, -1.0, 0.0]]
p0 = KMCProcess(coordinates=coordinates_p0,
                elements_before=['P','E'],
                elements_after=['E','P'],
                move_vectors=None,
                basis_sites=[0],
                rate_constant=1.)

# Setup a diffusion process to the right.
#                  center site      second site on 1
coordinates_p1 = [[0.0, 0.0, 0.0],[0.0, 1.0, 0.0]]
p1 = KMCProcess(coordinates=coordinates_p1,
                elements_before=['P','E'],
                elements_after=['E','P'],
                move_vectors=None,
                basis_sites=[0],
                rate_constant=1.)
# Setup a diffusion process to the top.
coordinates_p2 = [[0.0, 0.0, 0.0],[1.0, 0.0, 0.0]]
p2 = KMCProcess(coordinates=coordinates_p2,
                elements_before=['P','E'],
                elements_after=['E','P'],
                move_vectors=None,
                basis_sites=[0],
                rate_constant=1.)
# Setup a diffusion process to the right.
coordinates_p3 = [[0.0, 0.0, 0.0],[-1.0, 0.0, 0.0]]
p3 = KMCProcess(coordinates=coordinates_p3,
                elements_before=['P','E'],
                elements_after=['E','P'],
                move_vectors=None,
                basis_sites=[0],
                rate_constant=1.)

processes = [p0,p1,p2,p3]