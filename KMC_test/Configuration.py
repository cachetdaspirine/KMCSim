import numpy as np
import math
import copy

def distance(xy0,xy1):
    return np.sqrt((xy0[0]-xy1[0])**2+(xy0[1]-xy1[1])**2+(xy0[2]-xy1[2])**2)
def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper
moves_v = np.array([[1.,0.,0.],[-1.,0.,0.],[0.,1.,0.],[0.,-1.,0.]])
distance_set = {0,2.,truncate(np.sqrt(2),3),truncate(np.sqrt(5),3)}#,3,truncate(2*2**0.5,3)}
class Configuration:
    """
    configuration class on which the moves are applied.
    -> the positions of the particles are stored in the referential
    of the i=0 one:
    the class contain the required checks :
    -> check_excluded : returns False if any particles violate
                        the excluded volume requirement
    -> check_distances : returns False if any distances violate
                         the bonds distances requirement
    """
    def __init__(self,position_vector=None,Nparticles = 3):
        self.Nparticles=Nparticles
        # initial configuration with Nparticles aligned separated by an empty site
        if type(position_vector) == np.ndarray:
            self.position_v = position_vector
        else:
            self.position_v = np.array([[0.,2.*i,0.] for i in range(self.Nparticles)])
        self.initial_position = copy.copy(self.position_v)
    def get_tuple_position(self):
        return tuple(map(tuple, self.position_v))
    def reset_position(self):
        self.position_v = copy.copy(self.initial_position)
    def check_excluded(self):
        # check if no particle is within the neighboring of another one
        for xy0 in self.position_v:
            for xy1 in self.position_v:
                if (xy0 != xy1).any():
                    if sum(abs(xy0-xy1))<2:
                        return False
        return True
    def check_distances(self):
        for i in range(self.position_v.shape[0]-1):
            #for xy1 in self.position_v:
            xy0 = self.position_v[i]
            xy1 = self.position_v[i+1]
            if not truncate(distance(xy0,xy1),3) in distance_set:
                return False
        return True
    def apply_move(self,i,v):
        """
        apply the moving vector v to the ith particle
        then translate the system to have the 0th
        particle at point (0,0)
        """
        self.position_v[i] += v
        self.position_v-= self.position_v[0]
    def print_configuration(self):
        to_print = [['0 ' for _ in range(2*self.Nparticles+2)] for _ in range(2*self.Nparticles+2)]
        #to_print.fill(" 0 ")
        for ij in self.position_v:
            to_print[ij[0]-self.position_v[1,0]+3][ij[1]-self.position_v[1,1]+3] = "1 "
            #to_print[ij[0]+3][ij[1]+3] = "1 "
            # First row
        print(f"  ", end='')
        for j in range(8):
            print(f"| {j+1} ", end='')
        print(f"| ", end='')
        print()
        print(f'{(8*4+4)*"-"}')
        A=''
        for i in range(len(to_print)):
            print(f"{chr(65+i)} ", end='')
            for j in range(len(to_print[0])):
                if to_print[i][j]=='0 ':
                    print(f"| {' '} ", end='')
                else:
                    print(f"| {'1'} ", end='')
            print(f"| ", end='')
            print()
            print(f'{(8*4+4)*"-"}')
                #print(to_print[j][i],end='')
        print(A)