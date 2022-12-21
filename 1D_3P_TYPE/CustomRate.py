import numpy as np
import math
from KMCLib import *
from KMCLib.PluginInterfaces.KMCRateCalculatorPlugin import KMCRateCalculatorPlugin



class CustomPlugin(KMCRateCalculatorPlugin):
    def rate(self,coords,types_before,types_after,rate_constant,process_number,global_coordinate):
        moved_ID = self.configuration.movedAtomIDs()
        #print('moved_ID = '+str(moved_ID))
        if moved_ID.__len__()==1:
            return rate_constant*np.exp(Vs[types_after[0]](global_coordinate))
            """if types_after[0] =='A':
                return rate_constant*np.exp(VA(global_coordinate))
            elif types_after[0]=='B':
                return rate_constant*np.exp(VB(global_coordinate))
            elif types_after[0]=='C':
                return rate_constant*np.exp(VC(global_coordinate))"""
        if moved_ID.__len__()==2:
            E_before = sum([Vs[types_before[i]](global_coordinate+coords[i]) for i in range(types_before.__len__())])
            E_after = sum([Vs[types_after[i]](global_coordinate+coords[i]) for i in range(types_before.__len__())])
            return rate_constant*np.exp(-E_after+E_before)
        return 1.
    
    def cutoff(self):
        """ Determines the cutoff for this custom model """
        return 0.0