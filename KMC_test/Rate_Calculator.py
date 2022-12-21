import numpy as np
import math
from KMCLib import *
from KMCLib.PluginInterfaces.KMCRateCalculatorPlugin import KMCRateCalculatorPlugin

def distance(xy0,xy1):
    #print(xy0)
    #print(xy1)
    return np.sqrt((xy0[0]-xy1[0])**2+(xy0[1]-xy1[1])**2+(xy0[2]-xy1[2])**2)
def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper
distance_set = {2.,truncate(np.sqrt(2),3),truncate(np.sqrt(5),3)}

class CustomPlugin(KMCRateCalculatorPlugin):
    def rate(self,coords,types_before,types_after,rate_constant,process_number,global_coordinate):
        moved_ID = self.configuration.movedAtomIDs()
        #print('moved_ID = '+str(moved_ID))
        if moved_ID.__len__()==2:
            moved_coordinates = self.configuration.atomIDCoordinates()[moved_ID[0]] # 1 si pas actualisé avant
            #print('coordinates ='+str(moved_coordinates))
            bound_IDs = self.configuration.bound_IDs[moved_ID[0]]
            #print('bound_IDs = '+str(bound_IDs))
            Neighboring_coordinates = [self.configuration.atomIDCoordinates()[IDs] for IDs in bound_IDs]
            #print('neighboring coordinates = '+str(Neighboring_coordinates))
            for xy in Neighboring_coordinates:
                if not distance(xy,moved_coordinates) in distance_set :
                    return 10**(-10)
        #print('coordinate ='+str(self.configuration.atomIDCoordinates()[self.configuration.movedAtomIDs()]))
        return 1.
    