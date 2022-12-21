import numpy as np

from KMCLib import *
from KMCLib.PluginInterfaces.KMCAnalysisPlugin import KMCAnalysisPlugin

# create a custom measurement on the fly
class SingleTrajectory(KMCAnalysisPlugin):
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
        if not self.__track_type in configuration.possibleTypes():
            raise Error("The track type of the SingleTrajectory calculator is not one of the valid types of the configuration.")
        # compute the number of particles that are being tracked
        number_of_particle_tracked = configuration.particlesPerType()[configuration.possibleTypes()[self.__track_type]]
        # initialize a multi-dimensional array in which the trajectories will be saved
        # trajectories[:,time] all positions at a certain time
        # trajectories[particle,:] time evolution of a single particle
        self.__trajectories = np.zeros((number_of_particle_tracked,
                                        self.__step_max//self.__save_step,3),
                                       dtype=float)
        # stores the IDs of all the atoms we will track
        self.__IDs =np.where(np.array(configuration.atomIDTypes())==self.__track_type)[0]
        # initialize the trajectory to their initial position
        print(self.__IDs)
        self.__trajectories[:,0] = configuration.atomIDCoordinates()[self.__IDs]
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
        self.__trajectories[:,(step-1)//self.__save_step] = configuration.atomIDCoordinates()[self.__IDs]
    def finalize(self):
        """
        Finalize the trajectory saving
        nothing to do
        """

    def results(self):
        """
        Query function for the result.
        :returns: The result as a N_particle x step_tot//save_step
                  numpy array. All the trajectories are given in lattice number
                  units
        """
        return self.__trajectories
