import numpy as np

def Convert_to_vmd(in_name,out_name,polymer_size):
    input_file = open(in_name,'r')
    output_file = open(out_name,'w')
    #output_file.write(str(polymer_size)+'\n')
    lattice = ''
    repetition=''
    for i in range(10):
        # extract the lattice mesh
        ligne = input_file.readline()
        if 'a:' in ligne or 'b:' in ligne or 'c:' in ligne:
            lattice += ligne[3:-2]+" "
        # extract the repetition
        if 'REPETITIONS' in ligne:
            repetition = [float(i) for i in filter(None,ligne[11:].split(" "))]
    lattice =  np.array([float(i) for i in filter(None,lattice.split(" "))])
    
    lattice[0:3] = lattice[0:3]*repetition[0]
    lattice[3:6] = lattice[3:6]*repetition[1]
    lattice[6:9] = lattice[6:9]*repetition[2]
    #lattice = np.array2string(lattice, precision=2, separator=' ')
    
    for ligne in input_file:
        if 'STEP 0' in ligne:
            continue
        if 'STEP' in ligne:
            #output_file.write('\n')
            continue
        if 'TIME' in ligne:
            timestep = float(list(filter(None,ligne.split(' ')))[-1])
            output_file.write(str(polymer_size)+'\n')
            output_file.write('Lattice=\"'+str(lattice)[1:-1]+"\" Properties=species:S:1:pos:R:3 Time="+str(timestep)+"\n")
            
            #output_file.write(ligne)
        if 'P' in ligne:
            output_file.write(" ".join(ligne.split())+'\n')
    output_file.close()