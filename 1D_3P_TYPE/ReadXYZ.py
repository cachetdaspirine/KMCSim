def read_xyz_file(file_name):
  # open the .xyz file for reading
  with open(file_name, 'r') as f:
    # initialize a list to store the time steps
    time_steps = []

    # skip the header lines
    while True:
      line = f.readline()
      if line.startswith("PERIODICITY"):
        f.readline()
        break

    # read the file line by line
    while True:
      # read the first line containing the "STEP i" string
      line = f.readline()
      if not line:
        break
      step = int(line.split()[1])

      # read the second line containing the number of atoms
      num_atoms = int(f.readline().strip())

      # read the third line containing the "TIME t" string
      time = float(f.readline().split()[1])

      # initialize a list to store the coordinates for this time step
      coordinates = []

      # read the remaining lines, which contain the atom symbol, coordinates, and index
      for _ in range(num_atoms):
        line = f.readline().strip().split()
        symbol = line[0]
        x = float(line[1])
        y = float(line[2])
        z = float(line[3])
        index = int(line[4])
        coordinates.append((symbol, x, y, z, index))

      # add the time step data to the list
      time_steps.append((step, time, coordinates))

  # get the last time step
  step, time, coordinates = time_steps[-1]

  # sort the coordinates by x position
  coordinates.sort(key=lambda x: x[1])

  # print the atom types for the last time step sorted by x position
  print(f"STEP {step} (TIME {time}):")
  atom_types = [symbol if symbol!='S' else '_' for symbol, _, _, _, _ in coordinates]
  print(atom_types)