import numpy as np
import math
import copy

from Configuration import *

# Compute all the possible position vector
print('Start generating the configuration of middle particle in a polymer')
print('##################################################################')
# set that stores a unique sample of configuration
config_set_tot = set()
new_config_set = set()
config = Configuration()
# initialize the set with the first configuration
new_config_set.add(config.get_tuple_position())
# if all the position in new_config_set are already in config_set_top
# it means that we have all the configuration have been explored
while (new_config_set-new_config_set.intersection(config_set_tot)).__len__()!=0:
    config_set_tot |= new_config_set
    #print((new_config_set-new_config_set.intersection(config_set_tot)).__len__())
    previous_new_config_set = copy.copy(new_config_set)
    new_config_set = set()
    for position in previous_new_config_set:
        config = Configuration(np.array(position))
        for v in moves_v:
            for p in range(config.Nparticles):
                config.reset_position()
                config.apply_move(p,v)
                #config.print_configuration()
                if not config.check_distances() or not config.check_excluded():
                    continue
                new_config_set.add(config.get_tuple_position())
#print('Total configuration of three particles : '+str(config_set_tot.__len__()))

def test_if_middle_particle_moved(position1,position2):
    """
    Given two position vector, this function test if it is possible
    to go from one to the other by just moving the middle particle.
    """
    if position1[0]!=position2[0] or position1[2]!=position2[2]:
        return False
    if abs(position1[1][0]-position2[1][0])+abs(position1[1][1]-position2[1][1])==1:
        #print(np.append(np.array(position1),np.array([position2[1]]),axis=0))
        #print(position1,position2)
        return True
    return False

Coordinates_middle = list()
for position1 in config_set_tot:
    #old_size = Coordinates_middle.__len__()
    for position2 in config_set_tot:
        if test_if_middle_particle_moved(position1,position2):
            Coordinates_middle.append(np.append(np.array(position1),np.array([position2[1]]),axis=0))
    #if Coordinates_middle.__len__()-old_size==0:
    #    config = Configuration(np.array(position1))
    #    config.print_configuration()
Coordinates_middle=np.array(Coordinates_middle)
print('Total middle moves : Coordinates_middle.shape = '+str(Coordinates_middle.shape))
print('                     ------------------------')
#print('Start generating the configuration of side particle in a polymer')
# set that stores a unique sample of configuration
config_set_tot2 = set()
new_config_set = set()
config = Configuration(Nparticles=2)
# initialize the set with the first configuration
new_config_set.add(config.get_tuple_position())
# if all the position in new_config_set are already in config_set_top
# it means that we have all the configuration have been explored
while (new_config_set-new_config_set.intersection(config_set_tot2)).__len__()!=0:
    config_set_tot2 |= new_config_set
    #print((new_config_set-new_config_set.intersection(config_set_tot)).__len__())
    previous_new_config_set = copy.copy(new_config_set)
    new_config_set = set()
    for position in previous_new_config_set:
        config = Configuration(np.array(position),Nparticles=2)
        for v in moves_v:
            for p in range(config.Nparticles):
                config.reset_position()
                config.apply_move(p,v)
                #config.print_configuration()
                if not config.check_distances() or not config.check_excluded():
                    continue
                new_config_set.add(config.get_tuple_position())
#print('total number of two particles move : '+str(config_set_tot2.__len__()))
def test_if_2_particle_moved(position1,position2):
    """
    Given two position vector, this function test if it is possible
    to go from one to the other by just moving the middle particle.
    """
    if abs(position1[1][0]-position2[1][0])+abs(position1[1][1]-position2[1][1])==1:
        #print(np.append(np.array(position1),np.array([position2[1]]),axis=0))
        return True
    return False
Coordinates_extrem = list()
for position1 in config_set_tot2:
    #old_size = Coordinates_middle.__len__()
    for position2 in config_set_tot2:
        if test_if_2_particle_moved(position1,position2):
            Coordinates_extrem.append(np.append(np.array(position1),np.array([position2[1]]),axis=0))
    #if Coordinates_middle.__len__()-old_size==0:
    #    config = Configuration(np.array(position1))
    #    config.print_configuration()
Coordinates_extrem = np.array(Coordinates_extrem)
print('total extremities moves : Coordinates_extrem.shape = '+str(Coordinates_extrem.shape))
print('                          ------------------------')