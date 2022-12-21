from KMCLib import *

def make_interactions(kbc0,kac0,Aab=1.,Aac=1.,Abc=1.,kdiff=1.,kab0=1.,ChemicalMove=True):
    # Setup a diffusion process to the right.
    coordinates_p0 = [[1.0, 0.0, 0.0],[0.0, 0.0, 0.0]]
    p0A = KMCProcess(coordinates=coordinates_p0,
                    elements_before=['S','A'],
                    elements_after=['A','S'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kdiff)
    p0B = KMCProcess(coordinates=coordinates_p0,
                    elements_before=['S','B'],
                    elements_after=['B','S'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kdiff)
    p0C = KMCProcess(coordinates=coordinates_p0,
                    elements_before=['S','C'],
                    elements_after=['C','S'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kdiff)
    p0AB = KMCProcess(coordinates=coordinates_p0,
                    elements_before=['B','A'],
                    elements_after=['A','B'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kdiff)
    p0AC = KMCProcess(coordinates=coordinates_p0,
                    elements_before=['C','A'],
                    elements_after=['A','C'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kdiff)
    #coordinates_p0b = [[1.0, 0.0, 0.0],[0.0, 0.0, 0.0]]                
    p0BC = KMCProcess(coordinates=coordinates_p0,
                    elements_before=['C','B'],
                    elements_after=['B','C'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kdiff)
    # Setup a diffusion process to the left.
    coordinates_p1 = [[0.0, 0.0, 0.0],[-1.0, 0.0, 0.0]]
    p1A = KMCProcess(coordinates=coordinates_p1,
                    elements_before=['A','S'],
                    elements_after=['S','A'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kdiff)
    p1B = KMCProcess(coordinates=coordinates_p1,
                    elements_before=['B','S'],
                    elements_after=['S','B'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kdiff)
    p1C = KMCProcess(coordinates=coordinates_p1,
                    elements_before=['C','S'],
                    elements_after=['S','C'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kdiff)
    p1AB = KMCProcess(coordinates=coordinates_p1,
                    elements_before=['A','B'],
                    elements_after=['B','A'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kdiff)
    p1AC = KMCProcess(coordinates=coordinates_p1,
                    elements_before=['A','C'],
                    elements_after=['C','A'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kdiff)
    p1BC = KMCProcess(coordinates=coordinates_p1,
                    elements_before=['B','C'],
                    elements_after=['C','B'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kdiff)
    # Setup the chemical transition : 
    coordinates_chem = [[0.,0.,0.]]
    pAB = KMCProcess(coordinates=coordinates_chem,
                    elements_before=['A'],
                    elements_after=['B'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kab0*Aab)
    pBA = KMCProcess(coordinates=coordinates_chem,
                    elements_before=['B'],
                    elements_after=['A'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kab0)
    pAC = KMCProcess(coordinates=coordinates_chem,
                    elements_before=['A'],
                    elements_after=['C'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kac0*Aac)
    pCA = KMCProcess(coordinates=coordinates_chem,
                    elements_before=['C'],
                    elements_after=['A'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kac0)
    pBC = KMCProcess(coordinates=coordinates_chem,
                    elements_before=['B'],
                    elements_after=['C'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kbc0*Abc)
    pCB = KMCProcess(coordinates=coordinates_chem,
                    elements_before=['C'],
                    elements_after=['B'],
                    move_vectors=None,
                    basis_sites=[0],
                    rate_constant=kbc0)

    # Create the interactions object.
    if ChemicalMove:
        return KMCInteractions(processes=[p0A,p1A,p0B,p1B,p0C,p1C,p0AB,p0AC,p0BC,p1AB,p1AC,p1BC,
                                                pAB,pBA,pAC,pCA,pBC,pCB],
                                    implicit_wildcards=True)
    else :
        return KMCInteractions(processes=[p0A,p1A,p0B,p1B,p0C,p1C,p0AB,p0AC,p0BC,p1AB,p1AC,p1BC],                                                
                                    implicit_wildcards=True)
#interactions = KMCInteractions(processes=[p0A,p1A,p0B,p1B,p0C,p1C,p0AB,p0AC,p0BC,p1AB,p1AC,p1BC,pAB,pBA,pAC,pCA,pBC,pCB],
#                               implicit_wildcards=True)                               