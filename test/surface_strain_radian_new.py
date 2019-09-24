from ase.build import fcc111
from ase.constraints import FixAtoms, FixedPlane
from ase.visualize import view
from ase.io import write
from ase.io.vasp import write_vasp
from ase import Atoms, Atom
import numpy as np
import os, shutil
import argparse
from math import pi, cos, sin, ceil


# get conditions of system
pars = argparse.ArgumentParser()
pars.add_argument('-ang', type=str, nargs=3, default=['0', '30', '5'], help='MIN, MAX, INCREMENTS of theta. x axis is the criteria. ex) 0 30 5')
pars.add_argument('-ep', type=str, nargs=3, default=['-0.08','0.08','0.02'], help='MIN, MAX, INCREMENTS of epsilon(strain). ep=0 is pristine. ex) -0.08 0.08 0.02')
pars.add_argument('-atom', type=str, default='Cu')
pars.add_argument('-a', type=float, default=3.633287)
pars.add_argument('-size', type=int, nargs=3, default=[1,1,9], help="supercell size. ex) 1 1 9")
pars.add_argument('-vac', type=str, default='12.0', help='vacuum size in Angstrom')
args = pars.parse_args()
theta, epsilon, atom_name, lattice_a, supercell_size, vac_size= args.ang, args.ep, args.atom, args.a, args.size, float(args.vac)

theta=list(map(lambda x:float(x),theta))
epsilon=list(map(lambda x:float(x),epsilon))


def angle_ab(a,b,rad=True,acute=True):
      '''calculate internal angle between a vector (a1,a2,a3) and b vector(b1,b2,b3)'''
      from math import pi, cos, acos, sqrt
      a1,a2,a3=a; b1,b2,b3=b
      ab=acos((a1*b1+a2*b2+a3*b3)/(sqrt(a1**2+a2**2+a3**2)*sqrt(b1**2+b2**2+b3**2)))
      if acute==True:
            if ab > pi/2:
                  ab=pi-ab
      if rad==True:
            return ab
      else:
            return ab/pi*180




slab = fcc111(atom_name,a=lattice_a, size=supercell_size)
slab.center(axis=2, vacuum=vac_size)

cell_slab = slab.get_cell()

#write pristine surface which don't get any strain.
write_vasp('clean.vasp', slab, label='pristine surface', direct=True, sort=False, symbol_count = False, long_format=False,vasp5=True,ignore_constraints=False)
 
R1=np.arange(theta[0], theta[1], theta[2])
R1=np.append(R1,theta[1])
R2=np.arange(epsilon[0], epsilon[1], epsilon[2])
R2=np.append(R2,epsilon[1])

x1, x2, y1, y2 = cell_slab[0,0], cell_slab[0,1], cell_slab[1,0], cell_slab[1,1]

for q in R1:
      # we will use radian
      q_r=q*pi/180
      # Defining strain vector (t1,t2,0)
      t1=cos(q_r)*x1-sin(q_r)*x2
      t2=sin(q_r)*x1+cos(q_r)*x2
      for e in R2:
            # when using np.arange function, the value '0' in array is not recognized as exact 0 (ex. 1.387e-17). So we add 1, and then compare with exact 1.
            if e+1==1:
                  pass
            else:
                  q_x = angle_ab((t1,t2,0),(x1,x2,0),rad=True,acute=True)
                  q_y = angle_ab((t1,t2,0),(y1,y2,0),rad=True,acute=True)
                  x1_t = (1+e*cos(q_x))*x1
                  x2_t = (1+e*cos(q_x))*x2
                  y1_t = (1+e*cos(q_y))*y1
                  y2_t = (1+e*cos(q_y))*y2

                  position = slab.get_scaled_positions()
                  slab.set_cell([(x1_t,x2_t,0),(y1_t,y2_t,0),(0,0,cell_slab[2,2])])

                  x = format(q,'.3f')
                  y = format(e,'.3f')
                  mask = [atom.tag > 3 and atom.tag < 7 for atom in slab]
                  slab.set_constraint(FixAtoms(mask=mask))

                  slab.set_scaled_positions(position)
                  
                  path = str(str(x)+'_'+str(y)+'.vasp')
                  write_vasp(path, slab, label='strained surface', direct=True, sort= False, symbol_count = False, long_format=False,vasp5=True,ignore_constraints=False)
