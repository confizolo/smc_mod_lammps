import numpy as np
from numba import jit

def load_lammps(dump,natoms,elements, timestep = False):
    nelements = len(elements)
    
    if timestep:
        outarray = np.empty((0, natoms, nelements+1))
    else:   
        outarray = np.empty((0, natoms, nelements))

    index, time, flag = dump.iterator(0)
    while (flag+1):
        temp = np.array(dump.vecs(time,*elements)).T
        temp = np.pad(temp, [(0, natoms - temp.shape[0]), (0, 0)], mode="constant")
        if timestep:
            temp = np.pad(temp, [(0, 0), (1, 0)], mode="constant", constant_values=time)
        temp = temp.reshape(1,natoms,-1)
        outarray = np.vstack([outarray,temp])
        index, time, flag = dump.iterator(1)
    return outarray
