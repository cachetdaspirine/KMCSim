import numpy as np

from KMCLib import *
from KMCLib.PluginInterfaces.KMCAnalysisPlugin import KMCAnalysisPlugin

# create a custom measurement on the fly
class Density_thermal_av(KMCAnalysisPlugin):
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
        self.__track_type = track_type # type of the particles tracked has to be iterable
        self.__save_step = save_step # steps of measurement
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
        for type in self.__track_type:
            if not type in configuration.possibleTypes():
                raise TypeError("The track type of the SingleTrajectory calculator is not one of the valid types of the configuration.")
        # compute the number of particles that are being tracked
        #number_of_particle_tracked = configuration.particlesPerType()[configuration.possibleTypes()[self.__track_type]]
        # Initialize an array for each particle tracked
        self.__Density = np.zeros((self.__track_type.__len__(),configuration.lattice().repetitions()[0]),dtype=float)
        self.store_previous_state = np.zeros((self.__track_type.__len__(),configuration.lattice().repetitions()[0]),dtype=float)
        self.time_before = 0.
        
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
        # config.atomIDTypes() returns a tuple where each entry is a type.
        # thus np.argwhere(np.asarray(config.atomIDTypes())=='A')[:,0] returns the ID of every A atom
        # config.atomIDCoordinates() is an array of coordinate per ID we thus got all the coordinates of the A atoms
        # the last [:,0] means we only access the 0th coordinate
        # we then convert it to integer (because it's a float) which gives us the lattice site
        #print('beginning of the analysis, print the state before the move that has just been applied:')
        #print(self.store_previous_state)
        #print('Now output the state that has just been applied :')
        self.__Density+=self.store_previous_state*(time-self.time_before)
        self.store_previous_state.fill(0.)
        self.time_before = time
        for n,type in enumerate(self.__track_type):
            try:
                #print('n = ',n)
                #print(np.asarray(configuration.atomIDCoordinates()[np.argwhere(np.asarray(configuration.atomIDTypes())==type)[:,0]][:,0],dtype=int)
                self.store_previous_state[n,np.asarray(configuration.atomIDCoordinates()[np.argwhere(np.asarray(configuration.atomIDTypes())==type)[:,0]][:,0],dtype=int)]=1.
                #print(time-self.time_before)

            except IndexError:
                print(configuration.atomIDCoordinates()[np.argwhere(np.asarray(configuration.atomIDTypes())==type)[:,0]][:,0])
                raise
        
        #self.__Density[1,int(configuration.atomIDCoordinates()[np.argwhere(np.asarray(configuration.atomIDTypes())=='B')[:,0]][:,0])]+=1
        #self.__Density[2,int(configuration.atomIDCoordinates()[np.argwhere(np.asarray(configuration.atomIDTypes())=='C')[:,0]][:,0])]+=1
    def finalize(self):
        """
        Finalize the trajectory saving
        normalize the self.Density vector
        """
        self.__Density = self.__Density/self.time_before

    def results(self):
        """
        Query function for the result.
        :returns: The result as a N_particle x step_tot//save_step
                  numpy array. All the trajectories are given in lattice number
                  units
        """
        return self.__Density
