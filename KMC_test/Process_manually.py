from KMCLib import *
# Setup a diffusion process up-down, middle
#  *,  E,  *  ->   *,  P1, *
#  P0, P1, P2 ->   P0, E,  P2
#  *,  *,  *  ->   *,  *,  *
#                   P0             P1            P2           Moved
coordinates = [[0.0, 0.0, 0.0],[1.0, 0.0, 0.0],[2.0,0.0,0.0],[1.0,1.0,0.0]]
p0 = KMCProcess(coordinates=coordinates,
                elements_before=['P0','P1','P2','E'],
                elements_after=['P0','E','P2','P1'],
                move_vectors=None,
                basis_sites=[0],
                rate_constant=1.)
#Reverse
#  *,  E,  *  <-   *,  P1, *
#  P0, P1, P2 <-   P0, E,  P2
#  *,  *,  *  <-  *,  *,  *
#                   P0             P1            P2           Moved
coordinates = [[0.0, 0.0, 0.0],[1.0, 0.0, 0.0],[2.0,0.0,0.0],[1.0,1.0,0.0]]
p0R = KMCProcess(coordinates=coordinates,
                elements_before=['P0','E','P2','P1'],
                elements_after=['P0','P1','P2','E'],
                move_vectors=None,
                basis_sites=[0],
                rate_constant=1.)
# Setup a diffusion process up-down, middle
#  *,  E,  *  ->   *,  P1, *
#  P0, P1, P2 ->   P0, E,  P2
#  *,  *,  *  ->   *,  *,  *
#                   P0             P1            P2           Moved
coordinates = [[0.0, 0.0, 0.0],[1.0, 0.0, 0.0],[2.0,0.0,0.0],[1.0,-1.0,0.0]]
p1 = KMCProcess(coordinates=coordinates,
                elements_before=['P0','P1','P2','E'],
                elements_after=['P0','E','P2','P1'],
                move_vectors=None,
                basis_sites=[0],
                rate_constant=1.)
#Reverse
#  *,  E,  *  <-   *,  P1, *
#  P0, P1, P2 <-   P0, E,  P2
#  *,  *,  *  <-  *,  *,  *
#                   P0             P1            P2           Moved
coordinates = [[0.0, 0.0, 0.0],[1.0, 0.0, 0.0],[2.0,0.0,0.0],[1.0,-1.0,0.0]]
p1R = KMCProcess(coordinates=coordinates,
                elements_before=['P0','E','P2','P1'],
                elements_after=['P0','P1','P2','E'],
                move_vectors=None,
                basis_sites=[0],
                rate_constant=1.)
# Setup a diffusion process left-right, middle
#  *,  E,  *  ->   *,  P1, *
#  P0, P1, P2 ->   P0, E,  P2
#  *,  *,  *  ->   *,  *,  *
#                   P0             P1            P2           Moved
coordinates = [[0.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0,2.0,0.0],[1.0,1.0,0.0]]
p2 = KMCProcess(coordinates=coordinates,
                elements_before=['P0','P1','P2','E'],
                elements_after=['P0','E','P2','P1'],
                move_vectors=None,
                basis_sites=[0],
                rate_constant=1.)
#Reverse
#  *,  E,  *  <-   *,  P1, *
#  P0, P1, P2 <-   P0, E,  P2
#  *,  *,  *  <-  *,  *,  *
#                   P0             P1            P2           Moved
coordinates = [[0.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0,2.0,0.0],[1.0,1.0,0.0]]
p2R = KMCProcess(coordinates=coordinates,
                elements_before=['P0','E','P2','P1'],
                elements_after=['P0','P1','P2','E'],
                move_vectors=None,
                basis_sites=[0],
                rate_constant=1.)
# Setup a diffusion process left-right, middle
#  *,  E,  *  ->   *,  P1, *
#  P0, P1, P2 ->   P0, E,  P2
#  *,  *,  *  ->   *,  *,  *
#                   P0             P1            P2           Moved
coordinates = [[0.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0,2.0,0.0],[-1.0,1.0,0.0]]
p2 = KMCProcess(coordinates=coordinates,
                elements_before=['P0','P1','P2','E'],
                elements_after=['P0','E','P2','P1'],
                move_vectors=None,
                basis_sites=[0],
                rate_constant=1.)
#Reverse
#  *,  E,  *  <-   *,  P1, *
#  P0, P1, P2 <-   P0, E,  P2
#  *,  *,  *  <-  *,  *,  *
#                   P0             P1            P2           Moved
coordinates = [[0.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0,2.0,0.0],[-1.0,1.0,0.0]]
p2R = KMCProcess(coordinates=coordinates,
                elements_before=['P0','E','P2','P1'],
                elements_after=['P0','P1','P2','E'],
                move_vectors=None,
                basis_sites=[0],
                rate_constant=1.)


############################################################################
# Diffusion of the extremities
coordinates = [[0.0, 0.0, 0.0],[0.0, 1.0, 0.0],[1.0,2.0,0.0]]
p2R = KMCProcess(coordinates=coordinates,
                elements_before=['P0','E','P2','P1'],
                elements_after=['P0','P1','P2','E'],
                move_vectors=None,
                basis_sites=[0],
                rate_constant=1.)
