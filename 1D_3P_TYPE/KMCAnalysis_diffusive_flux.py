import numpy as np

from KMCLib import *
from KMCLib.PluginInterfaces.KMCAnalysisPlugin import KMCAnalysisPlugin

def I(x,lattice):
        """
        x is a 3D vector point
        This function works together with get_space_point. given a x,
        it returns the corresponding index of get_space_point.
        !!!! Only works for 1 basis point, and squared lattice!!!!
        """
        x = x- lattice.basis()[0] # remove the coordinate of the first lattice basis point
        Ix = x[0]/lattice.unitCell().cellVectors()[0][0]
        Iy = x[1]/lattice.unitCell().cellVectors()[1][1]
        Iz = x[2]/lattice.unitCell().cellVectors()[2][2]
        return int(Ix + Iy*lattice.repetitions()[0]+ Iz*lattice.repetitions()[0]*lattice.repetitions()[1])
# create a custom measurement on the fly
class DiffusiveFlux(KMCAnalysisPlugin):
    """
    Class for performing on-the-fly single trajectory.
    The trajectory is stored as a lattice trajectory
    """
    def __init__ (self,
                  track_type,
                  save_step,
                  step_max
                 ):
        """
        Initialize the measurement of a single trajectory
        """
        try:
            iter(track_type)
        except TypeError:
            raise TypeError("the track_type must be iterable, if it is a single value write it as a list of a single element.")
        if type(track_type)==str:
            raise TypeError("the track_type must be iterable, if it is a single value write it as a list of a single element.")
        self.__track_type = track_type
        self.__save_step = save_step
        self.__step_max = step_max
    def setup(self, step, time, configuration):
        """
        Recieves the setup call from the before the MC loop.

        :param step: The step number of the simulation.
        :type step: int

        :param time: The time of the simulation.
        :type time: float

        :param configuration: The configuration of the simulation.
        :type configuration: KMCConfiguration
        """
        # Make sure the track type is one of the possible types.
        for track_type in self.__track_type:
            if not track_type in configuration.possibleTypes():
                raise TypeError("The track type of the SingleTrajectory calculator is not one of the valid types of the configuration.")
        self.flux = {type : np.zeros(configuration.lattice().repetitions()[0],dtype=float) for type in self.__track_type}
        #self.N = {type : np.zeros(configuration.lattice().repetitions()[0],dtype=float) for type in self.__track_type}
        #self.time_before = 0.
    def registerStep(self, step, time, configuration):
        """
        Recieves the step call from the MC loop.

        :param step: The step number of the simulation.
        :type step: int

        :param time: The time of the simulation.
        :type time: float

        :param configuration: The configuration of the simulation.
        :type configuration: KMCConfiguration
        """
        if step<self.__step_max//100:
            # 1/100th of the simulation is dedicated to equilibration
            return
        if configuration.movedAtomIDs().__len__() == 2:
            for i in range(2):
                type = configuration.atomIDTypes()[configuration.movedAtomIDs()[i]]
                #print(self.__track_type)
                if type in self.__track_type:
                    XIF = I(configuration.atomIDCoordinates()[configuration.movedAtomIDs()[i]],configuration.lattice())        
                    XJF = I(configuration.atomIDCoordinates()[configuration.movedAtomIDs()[(i+1)%2]],configuration.lattice())                    
                    #print(time-self.time_before)
                    self.flux[type][XIF]+= 1#(XIF-XJF)#/(time-self.time_before)
                    self.flux[type][XJF]-= 1#XIF-XJF
                    #self.N[type][flux_index]+= 1
        #else:
        #    print(configuration.movedAtomIDs().__len__())
        self.time = time
    def finalize(self):
        """
        Finalize the trajectory saving
        nothing to do
        """
        for key in self.__track_type:
            self.flux[key] = self.flux[key]/(self.time*0.99) #we only used 99% of the moves
        #    for i in range(self.flux[key].shape[0]):
        #        if self.N[key][i]!=0:
        #            self.flux[key][i] = self.flux[key][i] / self.N[key][i]

    def results(self):
        """
        Query function for the result.
        :returns: The result as a N_particle x step_tot//save_step
                  numpy array. All the trajectories are given in lattice number
                  units
        """
        return self.flux