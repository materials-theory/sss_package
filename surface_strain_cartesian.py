from ase.build import fcc111
from ase.constraints import FixAtoms, FixedPlane
from ase.visualize import view
from ase.io import write
from ase.io.vasp import write_vasp
from ase import Atoms, Atom
import numpy as np
import os, shutil

slab = fcc111('Cu',a=3.633287, size=(1, 1, 9))
slab.center(axis=2, vacuum=12.0)

cell_slab = slab.get_cell()

write_vasp('clean', slab, label='pristine surface', direct=True, sort=False, symbol_count = False, long_format=False,vasp5=True,ignore_constraints=False)

for i in np.arange(0.970, 1.030, 0.005):
        for j in np.arange(0.970, 1.030, 0.005):
            x1 = i*cell_slab[0,0]
            x2 = i*cell_slab[0,1]

            y1 = j*cell_slab[1,0]
            y2 = j*cell_slab[1,1]

            position = slab.get_scaled_positions()
            slab.set_cell([(x1,x2,0),(y1,y2,0),(0,0,cell_slab[2,2])])

            x = format(i,'.3f')
            y = format(j,'.3f')
            mask = [atom.tag > 3 and atom.tag < 7 for atom in slab]
            slab.set_constraint(FixAtoms(mask=mask))

            slab.set_scaled_positions(position)
            
            path = str(str(x)+'_'+str(y)+'.vasp')
            write_vasp(path, slab, label='strained surface', direct=True, sort= False, symbol_count = False, long_format=False,vasp5=True,ignore_constraints=False)