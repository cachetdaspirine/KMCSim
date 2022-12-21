import numpy as np
from KMCLib import *
from KMCLib.PluginInterfaces.KMCRateCalculatorPlugin import KMCRateCalculatorPlugin
from KMCAnalysis_density_distribution import *
from KMCAnalysis_single_trajectory import *
from KMCAnalysis_diffusive_flux import *
import Process
import time
#from KMCAnalysis_single_trajectory import *

class OneDSimulation:
    def __init__(self,size):
        # Define a squared unit cell of size 1x1x1.
        cell_vectors = [[   1.000000e+00,   0.000000e+00,   0.000000e+00],
                        [   0.000000e+00,   1.000000e+00,   0.000000e+00],
                        [   0.000000e+00,   0.000000e+00,   1.000000e+00]]
        # Define a unique point within the unit cell
        basis_points = [[   0.500000e+00,   0.500000e+00,   0.500000e+00]]
        self.size=size
        try:
            unit_cell = KMCUnitCell(cell_vectors=cell_vectors,
                                    basis_points=basis_points)
        except:
            raise
        try:
            self.lattice = KMCLattice(unit_cell=unit_cell,
                                repetitions=(size,1,1),
                                periodic=(False, False, False))
        except:
            raise

    def get_space_point(self):
        """
        return the real space position of the basis points.
        This function use the unit cell and the repetition
        informations to generate a set of basis points throughout space.
        """
        res = list()
        for b in self.lattice.basis():
            res.append([b[0],b[1],b[2]])
            for d in range(3):
                for l in range(1,self.lattice.repetitions()[d]):
                    line = list()
                    for db in range(3):
                        line.append(b[db]+l*self.lattice.unitCell().cellVectors()[d][db])
                    res.append(line)
        return np.array(res)
    
    def I(self,x):
        """
        x is a 3D vector point
        This function works together with get_space_point. given a x,
        it returns the corresponding index of get_space_point.
        !!!! Only works for 1 basis point, and squared lattice!!!!
        """
        x = x- self.lattice.basis()[0] # remove the coordinate of the first lattice basis point
        Ix = x[0]/self.lattice.unitCell().cellVectors()[0][0]
        Iy = x[1]/self.lattice.unitCell().cellVectors()[1][1]
        Iz = x[2]/self.lattice.unitCell().cellVectors()[2][2]
        return int(Ix + Iy*self.lattice.repetitions()[0]+ Iz*self.lattice.repetitions()[0]*self.lattice.repetitions()[1])

    def set_up_interactions(self,VA,VB,VC,VAB,VAC,VBC,kac0,kbc0,Aab,Aac,Abc):
        Vs = {'A':VA,'B':VB,'C':VC,'S':lambda x:0.,
                    'AB':VAB,'BA':VAB,'AC':VAC,'CA':VAC,'BC':VBC,'CB':VBC}
        outer_class = self
        class CustomPlugin(KMCRateCalculatorPlugin):
            def initialize(self):
                self.outer_class = outer_class
            def rate(self,coords,types_before,types_after,rate_constant,process_number,global_coordinate):
                
                # chemical move:
                if types_before.__len__()==1:
                    return rate_constant* np.exp(-(Vs[types_after[0]+types_before[0]](global_coordinate+coords[0])-Vs[types_before[0]](global_coordinate+coords[0])))
                else:
                    E_before = sum([Vs[types_before[i]](global_coordinate+coords[i]) for i in range(types_before.__len__())])
                    E_after = sum([Vs[types_after[i]](global_coordinate+coords[i]) for i in range(types_before.__len__())])
                    return rate_constant*np.exp((-E_after+E_before)/2)            
                
            def cutoff(self):
                """ Determines the cutoff for this custom model """
                return 0.0
        try:
            self.interactions  = Process.make_interactions(kac0,kbc0,Aab,Aac,Abc)
            self.interactions.setRateCalculator(rate_calculator=CustomPlugin)
        except:
            raise

    def set_up_config_rho(self,rho):
        Nparticles=0
        k=0
        while Nparticles<1:
            if k>10**3:
                raise ValueError("cannot find a correct number of particles check rho value!")
            Nparticles = np.random.poisson(rho*self.size)# total number of particles in the system
            k+=1
        self.set_up_config(Nparticles)

    def set_up_config(self,Nparticles):
        """
        This function create an initial state with Nparticles
        in a 1D array of size self.size. The particles are spread
        randomly. All the particles are initially A's.
        """
        types = ['S']*(self.size)# solvent everywhere
        if Nparticles:
            self.Nparticles=Nparticles
        if self.Nparticles>self.size:
            raise ValueError("too many particles for the size of the system")
        for i in np.random.randint(0,self.size,size=self.Nparticles):
            types[i] = 'A'
        try:
            self.config = KMCConfiguration(lattice=self.lattice,
                                types=types,
                                possible_types=['S','A','B','C'])
        except:
            raise
    
    def set_up_specific_config(self,types,Nparticles):
        self.Nparticles=Nparticles
        try:
            self.config = KMCConfiguration(lattice=self.lattice,
                                types=types,
                                possible_types=['S','A','B','C'])
        except:
            raise

    def set_up_model(self):
        if hasattr(self,'config') and hasattr(self,'interactions'):
            self.model = KMCLatticeModel(configuration=self.config,
                                        interactions=self.interactions)
        else:
            raise ValueError('the config or the interaction is not defined. Please run first : set_up_config, and set_up_interactions.')
    
    def RUN(self,step_tot,dump_interval,analysis_interval,Single_Trajectory = False,Density=True,Flux=True):
        try:
            self.control_parameters = KMCControlParameters(number_of_steps=step_tot,
                                            dump_interval=dump_interval,
                                            analysis_interval=analysis_interval,
                                            seed=None)
        except:
            raise
        self.MyAnalysis = list()
        self.Name_to_index = dict() # a dictionnary that associate a name of analysis to the index in the above list
        n=0
        if Density:
            if self.control_parameters.analysisInterval()!=1:
                raise ValueError('the analysis interval has to be equal to 1')
            self.MyAnalysis.append(Density_thermal_av(['A','B','C'],
                    self.control_parameters.analysisInterval(),
                    self.control_parameters.numberOfSteps()))
            self.Name_to_index['Density'] = n
            n+=1
        if Single_Trajectory:
            self.MyAnalysis.append(SingleTrajectory('A',
                    self.control_parameters.analysisInterval(),
                    self.control_parameters.numberOfSteps()))
            self.Name_to_index['SingleTrajectory'] = n
            n+=1
        if Flux :
            self.MyAnalysis.append(DiffusiveFlux(['A','B','C'],
                                    self.control_parameters.analysisInterval(),
                                    self.control_parameters.numberOfSteps()))
            self.Name_to_index['DiffusiveFlux']=n
            n+=1
        
        t0 = time.perf_counter()
        try:
            self.model.run(self.control_parameters,trajectory_filename='test.py',
                        trajectory_type='xyz',
                        analysis = self.MyAnalysis)
        except:
            raise
        print('execution time : '+str(time.perf_counter()-t0))
    
    def get_density_profile(self):
        return self.get_space_point()[:,0],self.MyAnalysis[0].results()/self.Nparticles
    
    def get_diffusive_flux(self):
        return self.get_space_point()[:,0],self.MyAnalysis[self.Name_to_index['DiffusiveFlux']].results()
        
