from KMCLib import *

def make_interactions(kbc0,kac0,Aab=1.,Aac=1.,Abc=1.,kdiff=1.,kab0=1.,ChemicalMove=True):
    # Setup a diffusion process to the right.
    coordinates_p0 = [[0.0, 0.0, 0.0],[-1.0, 0.0, 0.0]]
    p0A = KMCProcess(coordinates=coordinates_p0,
                    elements_before=['A','S'],
                    elements_after=['S','A'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kdiff)
    p1A = KMCProcess(coordinates=coordinates_p0,
                        elements_before=['S','A'],
                        elements_after=['A','S'],
                        move_vectors=None,
                        basis_sites=[0],
                        rate_constant=kdiff)
    return KMCInteractions(processes=[p0A,p1A],implicit_wildcards=False)
