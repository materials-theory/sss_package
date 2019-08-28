from ase.build import fcc111
from ase.constraints import FixAtoms, FixedPlane
from ase.visualize import view
from ase.io import write
from ase.io.vasp import write_vasp
from ase import Atoms, Atom
import numpy as np
import os, shutil
import argparse
from math import pi, cos, acos, sqrt


# get conditions of system
pars = argparse.ArgumentParser()
pars.add_argument('-ang', type=str, nargs=3, default=['0', '30', '5'], help='MIN, MAX, INCREMENTS of theta. x axis is the criteria. ex) 0 30 5')
pars.add_argument('-ep', type=str, nargs=3, default=['-0.08','0.08','0.02'], help='MIN, MAX, INCREMENTS of epsilon(strain). ex) -0.08 0.08 0.02')
pars.add_argument('-atom', type=str, default='Cu')
pars.add_argument('-a', type=float, default=3.633287)
pars.add_argument('-size', type=int, nargs=3, default=[1,1,9], help="supercell size. ex) 1 1 9")
pars.add_argument('-vac', type=str, default='12.0', help='vacuum size in Angstrom')
args = pars.parse_args()
theta, epsilon, atom_name, lattice_a, supercell_size, vac_size= args.ang, args.ep, args.atom, args.a, args.size, float(args.vac)

theta=list(map(lambda x:float(x),theta))
epsilon=list(map(lambda x:float(x),epsilon))


slab = fcc111(atom_name,a=lattice_a, size=supercell_size)
slab.center(axis=2, vacuum=vac_size)

cell_slab = slab.get_cell()

#write pristine surface which don't get any strain.
write_vasp('clean.vasp', slab, label='pristine surface', direct=True, sort=False, symbol_count = False, long_format=False,vasp5=True,ignore_constraints=False)

R1=np.arange(theta[0], theta[1], theta[2])
R1=np.append(R1,theta[1])
R2=np.arange(epsilon[0], epsilon[1], epsilon[2])
R2=np.append(R2,epsilon[1])

for q in R1:
      for e in R2:
            # when using np.arange function, the value '0' in array is not recognized as exact 0 (ex. 1.387e-17). So we add 1, and then compare with exact 1.
            if e+1==1:
                  pass
            else:
                  # when we regard z=0 in both a, b vectors. We can easily get angle of instersection, angle_ab(rad)
                  angle_ab=acos((cell_slab[0,0]*cell_slab[1,0]+cell_slab[0,1]*cell_slab[1,1])/(sqrt(cell_slab[0,0]**2+cell_slab[0,1]**2)*sqrt(cell_slab[1,0]**2+cell_slab[1,1]**2)))
                  x1 = (1+e*cos(q*pi/180))*cell_slab[0,0]
                  x2 = (1+e*cos(q*pi/180))*cell_slab[0,1]
                  y1 = (1+e*cos(angle_ab-q*pi/180))*cell_slab[1,0]
                  y2 = (1+e*cos(angle_ab-q*pi/180))*cell_slab[1,1]

                  position = slab.get_scaled_positions()
                  slab.set_cell([(x1,x2,0),(y1,y2,0),(0,0,cell_slab[2,2])])

                  x = format(q,'.3f')
                  y = format(e,'.3f')
                  mask = [atom.tag > 3 and atom.tag < 7 for atom in slab]
                  slab.set_constraint(FixAtoms(mask=mask))

                  slab.set_scaled_positions(position)
                  
                  path = str(str(x)+'_'+str(y)+'.vasp')
                  write_vasp(path, slab, label='strained surface', direct=True, sort= False, symbol_count = False, long_format=False,vasp5=True,ignore_constraints=False)
