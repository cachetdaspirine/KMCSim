import numpy as np
import Coordinates
from KMCLib import *

def Generate_Processes(polymer_length):
    processes = list()
    for i in range(polymer_length-2):
        P0,P1,P2 = 'P'+str(i),'P'+str(i+1),'P'+str(i+2)
        for coordinates in Coordinates.Coordinates_middle :
            processes.append(KMCProcess(coordinates=coordinates,
                                        elements_before=[P0,P1,P2,'E'],
                                        elements_after=[P0,'E',P2,P1],
                                        move_vectors=None,
                                        basis_sites=[0],
                                        rate_constant=1.))
            processes.append(KMCProcess(coordinates=coordinates,
                                        elements_before=[P0,'E',P2,P1],
                                        elements_after=[P0,P1,P2,'E'],
                                        move_vectors=None,
                                        basis_sites=[0],
                                        rate_constant=1.))
    for coordinates in Coordinates.Coordinates_extrem :
        processes.append(KMCProcess(coordinates=coordinates,
                                    elements_before=['P0','P1','E'],
                                    elements_after=['E','P1','P0'],
                                    move_vectors=None,
                                    basis_sites=[0],
                                    rate_constant=1.))
        processes.append(KMCProcess(coordinates=coordinates,
                                    elements_before=['E','P0','P1'],
                                    elements_after=['P0','P1','E'],
                                    move_vectors=None,
                                    basis_sites=[0],
                                    rate_constant=1.))
        P_1 = 'P'+str(polymer_length-1)
        processes.append(KMCProcess(coordinates=coordinates,
                                    elements_before=[P_1,'P1','E'],
                                    elements_after=['E','P1',P_1],
                                    move_vectors=None,
                                    basis_sites=[0],
                                    rate_constant=1.))
        processes.append(KMCProcess(coordinates=coordinates,
                                    elements_before=['E',P_1,'P1'],
                                    elements_after=[P_1,'P1','E'],
                                    move_vectors=None,
                                    basis_sites=[0],
                                    rate_constant=1.))
    return processes