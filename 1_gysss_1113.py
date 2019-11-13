import subprocess
from ase.io import read, write
import os
import argparse
import numpy as np

##########################################################################
## -----------------------------About code----------------------------- ##
## Giyeok Lee, Ethan                                                    ##
## sss_package of MTG (https://github.com/materials-theory/sss_package) ##
## Version : 2019.11.05 (By Giyeok Lee)                                 ##
## Revised : 2019.11.05 (By Giyeok Lee)                                 ##
## -------------------------------------------------------------------- ##
##########################################################################


############################################################################################
## -----------------------------In 'input' directory------------------------------------- ##
## INCAR1 : VC-opt                                                                        ##
## INCAR2 : constrained VC-opt. (z-optimization)                                          ##
## INCAR3 : geo-opt for slab model.                                                       ##
## INCAR4 : scf calculation                                                               ##
## INCAR5 : non-scf for Workfunction Calculation.                                         ##
## + POSCAR, POTCAR, KPOINTS(or KSPACING tag), and other compulsory input files           ##
## -------------------------------------------------------------------------------------- ##
############################################################################################


# needed to be modified
# KOHN
z_opt="mpirun -np $SLURM_NTASKS /home/giyeok/1_program/VASP_constrained/vasp_std_z_fix > stdout.log"
norm_vasp="mpirun -np $SLURM_NTASKS /TGM/Apps/VASP/bin/5.4.4/NORMAL/vasp_5.4.4_GRP7_NORMAL_20170903.x > stdout.log"

#nurion1
# z_opt="mpirun sdf > stdout.log"
# norm_vasp="mpirun /home01/x1718a09/program/vasp.5.4.4/bin/vasp_std > stdout.log"


def read_POSCAR(file_name="POSCAR"):
	'''Read structure files of VASP program.
    [Input] : file's name. Default="POSCAR"
    [Output] : system_name, lattice, cell_vec, species, natoms, coord_type, coord, Selective, selective_tags
	|-> system_name : name of system. @str
    |-> lattice : scale @float
    |-> cell_vec: [x1, y1, z1], [x2, y2, z2], [x3, y3, z3]] @np.array(dtype='d')
    |-> species: [atomt1, atom2, ...] @list(str)
    |-> natoms: [number_of_atom1, number_of_atom2, ...] @list(int)
	|-> coord_type: "DIRECT" / "CARTESIAN" @str
	|-> coord : coordination of each atoms @np.array(dtype='d')
    |-> Selective : True / False @bool
    |-> selective_tags : [[T/F, T/F, T/F], [T/F, T/F, T/F], ...] @np.array(dtype=str)
    '''
	ip=open(file_name, 'r')
	system_name=ip.readline().strip()
	lattice=float(ip.readline().strip())
	#cell vector
	cell_vec=[]
	for i in range(3):
	    cell_vec.append(ip.readline().split())
	cell_vec=np.array(cell_vec,dtype="d")

	# atomic information (species, number)
	species=ip.readline().split()
	nspecies=len(species)
	natoms=list(map(lambda x:int(x),ip.readline().split()))
	tot_natoms=sum(natoms)

	# Selective Dynamics & Coord_Type
	what1=ip.readline().strip()
	if what1.upper().startswith("S"):
		Selective=True
		what2=ip.readline().strip()
		if what2.upper().startswith("D"):
			Cartesian=False
		elif what2.upper().startswith("F"):
			Cartesian=False
		else:
			Cartesian=True
	else:
		Selective=False
		if what1.upper().startswith("D"):
			Cartesian=False
		elif what1.upper().startswith("F"):
			Cartesian=False
		else:
			Cartesian=True
	if Cartesian:
		coord_type="CARTESIAN"
	else:
		coord_type="DIRECT"

	# Coordination
	coord=[]
	selective_tags=[]
	if not(Selective):
		for i in range(tot_natoms):
			coord.append(ip.readline().split())
		try:
			coord=np.array(coord,dtype="d")
		except:
			raise IOError("POSCAR file is weird. In [%s] file, Maybe number of coordinations are not matched with total number of atoms"%(file_name))
	else:
		for i in range(tot_natoms):
			line=ip.readline().split()
			coord.append(line[0:3])
			selective_tags.append(line[3:])
		try:
			coord=np.array(coord,dtype="d")
		except:
			raise IOError("POSCAR file is weird. In [%s] file, Maybe number of coordinations are not matched with total number of atoms"%(file_name))
		selective_tags=np.array(selective_tags,dtype=str)
	ip.close()
	return system_name, lattice, cell_vec, species, natoms, coord_type, coord, Selective, selective_tags


def write_POSCAR(output_file, system_name, lattice, cell_vec, species, natoms, coord_type, coord, Selective=False, selective_tags=[]):
	'''Write structure files of VASP program. Not including convert function here. dir2car and car2dir function have it.
    [Input] : output_file(name of output file), lattice, cell_vec, species, natoms, coord_type, coord, Selective, selective_tags
	|-> system_name : name of system. @str
    |-> lattice : scale @float
    |-> cell_vec: [x1, y1, z1], [x2, y2, z2], [x3, y3, z3]] @np.array(dtype='d')
    |-> species: [atomt1, atom2, ...] @list(str)
    |-> natoms: [number_of_atom1, number_of_atom2, ...] @list(int)
	|-> coord_type: "DIRECT" / "CARTESIAN" @str
	|-> coord : coordination of each atoms @np.array(dtype='d')
    |-> Selective : True / False @bool
    |-> selective_tags : [[T/F, T/F, T/F], [T/F, T/F, T/F], ...] @np.array(dltype=str)
    [Output] : VASP format structure file.
    '''
	op=open(output_file,'w')
	print(system_name,file=op)
	print("   %16.14f"%(lattice),file=op)
	for i in range(3):
	    print("    %19.16f   %19.16f   %19.16f"%(cell_vec[i][0],cell_vec[i][1],cell_vec[i][2]),file=op)

	print("",end=" ",file=op)
	for i in range(len(species)):
	    print("%4s"%(species[i]),end="",file=op)
	print("",file=op)

	for i in range(len(natoms)):
	    print("%6d"%(natoms[i]),end="",file=op)
	print("",file=op)

	if Selective:
	    print("Selective dynamics",file=op)
	if coord_type=="CARTESIAN":
	    print("Cartesian",file=op)
	else:
	    print("Direct",file=op)

	if Selective:
	    for (x1,y1) in zip(coord,selective_tags):
	        for i in x1:
	            print("%20.16f"%(i),end="",file=op)
	        if len(y1)==1:
	        	y1=y1[0].split()
	        for p in y1:
	            print("%4s"%(p),end="",file=op)
	        print("",file=op)
	else: 
	    for x1 in coord:
	        for i in x1:
	            print("%20.16f"%(i),end="",file=op)
	        print("",file=op)
	op.close()
	print("(+) Structure file [%s] was generated"%(output_file))

def dir2car(coord,cell_vec):
    '''Direct to Cartesian
    [Input] : coord, cell_vec
    |-> coord : coordination. @np.array(dtype="d")
    |-> cell_vec : cell vector. @np.array(dtype='d')
    [Output] : new coordination @np.array(dtype='d')
    '''
    return np.matmul(coord,cell_vec)

def car2dir(coord,cell_vec):
    '''Cartesian to Direct
    [Input] : coord, cell_vec
    |-> coord : coordination. @np.array(dtype='d')
    |-> cell_vec : cell vector. @np.array(dtype='d')
    [Output] : new coordination @np.array(dtype='d')
    '''
    inv_cell_vec=np.linalg.inv(cell_vec)
    return np.matmul(coord,inv_cell_vec)


def extract(object="OUTCAR",property="energy"):
	energy="cat OUTCAR | grep \"energy without entropy\" | tail -1 | awk '{printf \"%10.9f\", $5}'"
	fermi="grep fermi OUTCAR | tail -1 | awk '{printf\"%20.8f \",$3}'"
	volume="cat OUTCAR | grep volume | tail -1 | awk '{printf \"%20.8f\",$5}'"
	if property=="energy":
		# extract(property="energy")
		i=float(subprocess.check_output(energy,shell=True))
	elif property=="fermi":
		# extract(property="fermi")
		i=float(subprocess.check_output(fermi,shell=True))
	elif property=="volume":
		# extract(property="volume")
		i=float(subprocess.check_output(volume,shell=True))
	return i

def run_vasp(vasp_command):
	'''i=0 means exit 0, which is terminated normally. if not, it should be calculated again'''
	i=subprocess.call(vasp_command,shell=True)
	return i
#cal1_exit=run_vasp(z_opt)


def mkdir(dirname):
    directory = os.path.dirname(dirname+"/")
    if not os.path.exists(directory):
        os.makedirs(directory)


	#####

	  # def is_converge(self):
   #      if self.convergence_info is None or len(self.convergence_info) < 3:
   #          return False

   #      energies = [x['free_energy'] for x in self.convergence_info]
   #      if len(energies) > 2 and abs(max(energies[-3:]) - min(energies[-3:])) >= self.energy_tolerance:
   #          return False
   #      else:
   #          return True


def angle_ab(a,b,rad=True,acute=True):
	'''
	Calculate internal angle between first vector (a1,a2,a3) and second vector(b1,b2,b3)
	[Input] a, b, rad, acute
	|---> a : first vector. @list/tuple/np.array
	|---> b : second vector. @list/tuple/np.array
	|---> rad : True means output's unit will be radian. False means output's unit will be degree @bool.
	|---> acute : acute literally means the acute angle. if True : output doens't exceed 90 deg, False doesn't care. @bool
	[Output] : Angles, unit is degree(when rad==False) or radian(when rad==True)
	'''
	from math import pi, cos, sin, acos, sqrt
	a1,a2,a3=a; b1,b2,b3=b
	ab=acos((a1*b1+a2*b2+a3*b3)/(sqrt(a1**2+a2**2+a3**2)*sqrt(b1**2+b2**2+b3**2)))
	if acute==True:
	    if ab > pi/2:
	          ab=pi-ab
	if rad==True:
	    return ab
	else:
	    return ab/pi*180



def biaxial_strain(cell_vec,theta,epsilon,fractional=True):
	'''
	[Input] : cell vector, angle(theta), strain(epsilon)
	|---> theta : Angles between applied strain's axis and x axis. Units are degree. 0 means parallel to x axis @float
	|---> epsilon : Length of strain along strain axis. @float
	|---> fractional : False when the epsilon value is the exact length. True when epsilon value is ratio value (compare to length of diagonal line).
			|----> In both case 'fractional=True and fractional=False', epsilon=0 is pristine
	[Output] : new cell vector @np.array(dtype='d')
	'''
	from math import pi, cos, sin, acos, sqrt
	a1, b1, a2, b2 = cell_vec[0,0], cell_vec[0,1], cell_vec[1,0], cell_vec[1,1]

	def length(a):
		'''a : np.array, dtype=float'''
		temp=0
		for i in a:
			temp+=i**2
		return temp**(0.5)

	def cartesian_strain(cell_vec,xdiff,ydiff):
		'''
		[Input] : cell vector, difference along x and y axis.
		|---> cell_vec : cell vector @np.array(dtype='d')
		|---> xdiff : applied strain along x directions. 0.00 means no applied strain @float
		|---> ydiff : applied strain along y directions. 0.00 means no applied strain @float
		[Output] : new cell_vector @np.array(dtype='d')
		'''
		cell_vec[0] *= (xdiff+1)
		cell_vec[1] *= (ydiff+1)
		return cell_vec

	theta_r=theta*pi/180

	# Define New Diagonal vector, as N_Diag
	Diag=(cell_vec[0]+cell_vec[1])[:-1]
	l_Diag=length(Diag)

	if fractional==True:
		epsilon*=(l_Diag)

	N_Diag=Diag+np.array([epsilon*cos(theta_r),epsilon*sin(theta_r)])


	# Define new cell vectors.
	New_x1=(-b2/a2*N_Diag[0]+N_Diag[1])/(b1/a1-b2/a2)
	New_y1=b1/a1*New_x1
	New_x2=(-b1/a1*N_Diag[0]+N_Diag[1])/(b2/a2-b1/a1)
	New_y2=b2/a2*New_x2

	New_X=length(np.array([New_x1,New_y1]))
	New_Y=length(np.array([New_x2,New_y2]))


	Xdiff=(New_X-length(cell_vec[0]))/length(cell_vec[0])
	Ydiff=(New_Y-length(cell_vec[1]))/length(cell_vec[1])
	N_cell_vec=cartesian_strain(cell_vec,Xdiff,Ydiff)
	return N_cell_vec




#def biaxial_stress(cell_vec,theta,sigma):
#	'''
#	[Input] : cell vector, angle(theta), sigma(stress)
#	|---> theta : Angles between applied stress's axis and x axis. Units are degree. 0 means along x axis @float
#	|---> sigma : Amount of applied stress along stress axis. @float
#	[Output] : new cell vector @np.array(dtype='d')
#	'''
#	# from https://github.com/materials-theory/sss_package/blob/master/surface_strain_radian_v0924.py
#	from math import pi, cos, sin, acos, sqrt
#	x1, x2, y1, y2 = cell_vec[0,0], cell_vec[0,1], cell_vec[1,0], cell_vec[1,1]
#    # we will use radian
#	theta_r=theta*pi/180
#
#	# Defining stress vector (t1,t2,0)
#	t1=cos(theta_r)*x1-sin(theta_r)*x2
#	t2=sin(theta_r)*x1+cos(theta_r)*x2
#	q_x = angle_ab((t1,t2,0),(x1,x2,0),rad=True,acute=True)
#	q_y = angle_ab((t1,t2,0),(y1,y2,0),rad=True,acute=True)
#	x1_t = (1+stress*cos(q_x))*x1
#	x2_t = (1+stress*cos(q_x))*x2
#	y1_t = (1+stress*cos(q_y))*y1
#	y2_t = (1+stress*cos(q_y))*y2
#	cell_vec = np.array([(x1_t,x2_t,0),(y1_t,y2_t,0),(0,0,cell_vec[2,2])],dtype='d')
#	return cell_vec


def making_slab(input_file,dlayer,nlayer,flayer,vacuum_length):
    '''Making slabs along 001 direction.'''
    POSCAR=read_POSCAR(input_file)
    z=POSCAR[2][2][2]
    mlayer=(nlayer-flayer)//2
    def atomindex(index,natoms):
        '''To check the atom species by index. If atom species is ['I','Cu'], 'I' will be index 1 and 'Cu' will be index 2
        '''
        if index>sum(natoms):
            raise IOError("index exceeds the total number of atoms. please check about it")
        temp=[]
        for i in natoms:
            temp.append(i)
        for i in range(len(temp)):
            temp[i]+=sum(natoms[:i])
            if index<=temp[i]:
                return i+1



    if POSCAR[5]=="CARTESIAN":
        coord_car=POSCAR[6]
    else:
        coord_car=dir2car(POSCAR[6],POSCAR[2])

    ref=np.argsort(coord_car[:,2])


    n=1
    SB_coord=dict()
    SB_index=dict()


    temp=0
    for i in ref:
        if temp==0:
            SB_coord[n]=[coord_car[i]]
            SB_index[n]=np.array([i])
            temp+=1
        else:
            if coord_car[i][2]-np.mean(np.array(SB_coord[n])[:,2])<dlayer:
                SB_coord[n].append(coord_car[i])
                SB_index[n]=np.append(SB_index[n],i)
                if temp==len(ref)-1:
                    SB_coord[n]=np.array(SB_coord[n])
            else:
                SB_coord[n]=np.array(SB_coord[n])
                n+=1
                SB_coord[n]=[coord_car[i]]
                SB_index[n]=np.array([i])
                if temp==len(ref)-1:
                    SB_coord[n]=np.array(SB_coord[n])
            temp+=1

    # Interlayer spacing
    # [2nd-1st, 3rd-2nd, ... , 1st-topmost layer]
    IS_1=np.array([])
    for i in SB_coord:
        IS_1=np.append(IS_1,np.mean(SB_coord[i][:,2]))
    IS_2=IS_1[1:]
    IS_2=np.append(IS_2,IS_1[0]+z)
    IS=IS_2-IS_1

    if IS[-1]>np.mean(IS[:-2])+5:
        raise IOError('Reading POSCAR may be not a bulk. Please use bulk POSCAR.')


    Wcoord=dict()
    Wstags=dict()


    if n>=nlayer:
        for i in range(1,nlayer+1):
            for p in SB_index[i]:
                ats=atomindex(p+1,POSCAR[4])
                if ats not in Wcoord:
                    Wcoord[ats]=[]
                Wcoord[ats].append(coord_car[p])
                if i>mlayer and i<=mlayer+flayer:
                    if ats not in Wstags:
                        Wstags[ats]=[]
                    Wstags[ats].append(['F   F   F'])
                else:
                    if ats not in Wstags:
                        Wstags[ats]=[]
                    Wstags[ats].append(['T   T   T'])
    else:
        # when the number of needed layers is larger than original number of layers, This case will be dominant case.
        # n<nlayer
        for i in range(1,n+1):
            for p in SB_index[i]:
                ats=atomindex(p+1,POSCAR[4])
                if ats not in Wcoord:
                    Wcoord[ats]=[]
                Wcoord[ats].append(coord_car[p])
                if i>mlayer and i<=mlayer+flayer:
                    if ats not in Wstags:
                        Wstags[ats]=[]
                    Wstags[ats].append(['F   F   F'])
                else:
                    if ats not in Wstags:
                        Wstags[ats]=[]
                    Wstags[ats].append(['T   T   T'])

        for i in range(n+1,nlayer+1):
            quotient, remainder = i//n, i%n

            if remainder+1==1:
                if quotient+1==1:
                    quotient=0
                else:
                    quotient-=1
                remainder=n

            for p in SB_index[remainder]:
                ats=atomindex(p+1,POSCAR[4])
                if ats not in Wcoord:
                    Wcoord[ats]=[]
                temp=np.array([coord_car[p][0],coord_car[p][1]])
                temp=np.append(temp,coord_car[p][2]+z*quotient)
                # 여기서 coord_car[p]에 z방향으로 넣어주면 되는데, 원본이 바뀌지 않게 복사해서 넣어두자.
                Wcoord[ats].append(temp)            

                if i>mlayer and i<=mlayer+flayer:
                    if ats not in Wstags:
                        Wstags[ats]=[]
                    Wstags[ats].append(['F   F   F'])
                else:
                    if ats not in Wstags:
                        Wstags[ats]=[]
                    Wstags[ats].append(['T   T   T'])            




    Output_coord=[]
    Output_stags=[]
    Output_natoms=[]
    for i in range(1,len(Wcoord)+1):
        for p in Wcoord[i]:
            Output_coord.append(p)
        Output_natoms.append(len(Wcoord[i]))
    for i in range(1,len(Wstags)+1):
        for p in Wstags[i]:
            Output_stags.append(p)

    Output_coord=np.array(Output_coord)
    Output_natoms=np.array(Output_natoms)
    Output_stags=np.array(Output_stags)

    POSCAR[2][2][2]=max(Output_coord[:,2])+vacuum_length

    dir_coord=car2dir(Output_coord,POSCAR[2])
    write_POSCAR("POSCAR", POSCAR[0], POSCAR[1], POSCAR[2], POSCAR[3], Output_natoms, "DIRECT", dir_coord, Selective=True, selective_tags=Output_stags)


#################################### main part (prepare input files (INCAR1~5, KPOINTS, POSCAR, POTCAR, etc) in input/ directory.)

def report(cal, log, index):
	if cal==0:
		print("Calculation %1d finished."%(index), file=log)
	else:
		print("Calcuation %1d"%(index), file=log)

# (1) from optimized bulk (VC-opt), we got Energy <- INCAR1
def first_cal(log):

	# VASP Calculation
	index=1
	cal=run_vasp(norm_vasp)
	report(cal,log,index)

	E1_bulk=extract(property="energy")
	CONTCAR=read_POSCAR("CONTCAR")
	print("E1_bulk = [%5.5f]"%(E1_bulk),file=log)

	return cal, E1_bulk, CONTCAR


# (2) constrained bulk
def second_cal(in_POSCAR, log, theta, epsilon):
	# Applying biaxial strain.

	# Reading POSCAR
	POSCAR=read_POSCAR(in_POSCAR)
	if POSCAR[5]!="DIRECT":
		POSCAR[6]=car2dir(POSCAR[6],POSCAR[2])

	# Applying cartesian strain and make new POSCAR
	cell_vec=biaxial_strain(POSCAR[2], theta, epsilon)
	write_POSCAR("POSCAR", POSCAR[0], POSCAR[1], cell_vec, POSCAR[3], POSCAR[4], "DIRECT", POSCAR[6], POSCAR[7], POSCAR[8])

	# VASP Calculation
	index=2
	cal=run_vasp(z_opt)
	report(cal,log,index)

	E2_stbulk=extract(property="energy")
	print("E2_stbulk = [%5.5f]"%(E2_stbulk),file=log)

	CONTCAR=read_POSCAR("CONTCAR")
	return cal, E2_stbulk, CONTCAR

# (3) making slab from (2)
def third_cal(in_POSCAR, log, dlayer, nlayer, flayer, vacuum_length=10, slab=True):
	''' making slab
	[Slab] :  If slab==True, making slab using vacuum_length.
	|      :  If you want calculation from already-made slab POSCAR, use slab==False
	|---> Until now, Making asymmetric slab is only available.
	'''
	making_slab(in_POSCAR,dlayer,nlayer,flayer,vacuum_length)

	# VASP Calculation
	index=3
	cal=run_vasp(norm_vasp)
	report(cal,log,index)

	E3_slab=extract(property="energy")
	print("E3_slab = [%5.5f]"%(E3_slab),file=log)
	CONTCAR=read_POSCAR("CONTCAR")
	return cal, E3_slab, CONTCAR

# (4) scf calculation to make CHGCAR and WAVECAR
def forth_cal(log):

	# VASP Calculation
	index=4
	cal=run_vasp(norm_vasp)
	report(cal,log,index)

	E4_scf=extract(property="energy")
	print("E4_scf = [%5.5f]"%(E4_scf),file=log)
	return cal, E4_scf


# (5) non-scf calculation to make LOCPOT file to calculate Surface Work Function.
def fifth_cal(log):

	# VASP Calculation
	index=5
	cal=run_vasp(norm_vasp)
	report(cal,log,index)

	E5_nonscf=extract(property="energy")
	print("E5_nonsef = [%5.5f]"%(E5_nonscf), file=log)
	return cal, E5_nonscf



## Actual running parts.
pars = argparse.ArgumentParser()
pars.add_argument('-m', type = int, default=1, help='Choose which step do you want to start from. * For example, 1 (Default) : All (cal 1~5) / 2 : cal 2~5 ...')
pars.add_argument('-p', type = int, default=4, help='Choose the property you want to calculate. 1 : surface tension, 2 : workfunction, 3 : d-band center, 4 : All')
pars.add_argument('-vac', type = float, help="Choose the length of vacuum.", default=15.00)
pars.add_argument('-ang', type=str, nargs=3, default=['0', '0', '0'], help='MIN, MAX, INCREMENTS of theta. x axis is the criteria. ex) 0 30 5')
pars.add_argument('-ep', type=str, nargs=3, default=['0','0','0'], help='MIN, MAX, INCREMENTS of epsilon(strain). ep=0 is pristine. ex) -0.08 0.08 0.02')
pars.add_argument('-dlayer', type=float, default=0.5, help="Choose the width of layer to iterate.")
pars.add_argument('-nlayer', type=int, default=13, help="Choose the number of layers when making slabs.")
pars.add_argument('-flayer', type=int, default=3, help="Choose the number of center layers to fix (do not move)")
# pars.add_argument('-xst', type=str, nargs=3, default=['1','1','1'], help='MIN, MAX, INCREMENTS of applied strain along x axis. ex) 0.97 1.03 0.03')
# pars.add_argument('-yst', type=str, nargs=3, default=['1','1','1'], help='MIN, MAX, INCREMENTS of applied strain along y axis. ex) 0.97 1.03 0.03')
args = pars.parse_args()
mode, property_mode, vacuum_length, theta, epsilon, dlayer, nlayer, flayer = args.m, args.p, args.vac, args.ang, args.ep, args.dlayer, args.nlayer, args.flayer
# xstrain, ystrain=args.xst, args.yst

theta=list(map(lambda x:float(x),theta))
epsilon=list(map(lambda x:float(x),epsilon))

if theta[2]+1==1:
	theta=np.array([theta[0]])
else:
	theta_temp=np.arange(theta[0],theta[1],theta[2])
	theta=np.append(theta_temp,theta[1])
if epsilon[2]+1==1:
	epsilon=np.array([epsilon[0]])
else:
	epsilon_temp=np.arange(epsilon[0],epsilon[1],epsilon[2])
	epsilon=np.append(epsilon_temp, epsilon[1])

# xstrain=list(map(lambda x:float(x),xstrain))
# ystrain=list(map(lambda x:float(x),ystrain))
# xstrain_temp=np.arange(xstrain[0],xstrain[1],xstrain[2])
# xstrain=np.append(xstrain_temp,xstrain[1])
# ystrain_temp=np.arange(ystrain[0],ystrain[1],ystrain[2])
# ystrain=np.append(ystrain_temp,ystrain[1])


if not os.path.exists("./input/"):
	mkdir("process")
	log=open("process/[sss]calulation_log",'w')
	print("input files are not existed", file=log)
	log.close()
	raise IOError('''Can't find input files.
		in running directory,n

		directory/
			|__input/
				|__INCAR 1~5
				|__POSCAR
				|__KPOINTS
				|__POTCAR
				|__other input files (kernel, etc...)
			|__job_file (usually .sh format) <-- Which contains commands which executing this python scripts.

		this is right usage. Please make input files in 'input/' directory.
		''')
else:
	if mode==1:
		mkdir("1_bulk")
		os.chdir("1_bulk")
		subprocess.call("cp ../input/* ./; mv INCAR1 INCAR; rm INCAR2 INCAR3 INCAR4 INCAR5",shell=True)
		blog=open("[sss]calculation_log",'w')
		cal1, E1_bulk, CONTCAR_1 = first_cal(blog)
		blog.close()
		os.chdir("../")

	for th in theta:
		for ep in epsilon:

			if ep+1==1:
				process="process_clean"
				if os.path.exists(process):
					pass
				mkdir(process)
				log=open("%s/[sss]calulation_log"%(process),'w')
			else:
				process="process_%.2f_%.2f"%(th, ep)
				mkdir(process)
				log=open("%s/[sss]calulation_log"%(process),'w')

			if mode==1:
				mkdir("%s/2_constrained_bulk"%(process))
				os.chdir("%s/2_constrained_bulk"%(process))
				subprocess.call("cp ../../input/* ./; mv INCAR2 INCAR; rm INCAR1 INCAR3 INCAR4 INCAR5 POSCAR; cp ../../1_bulk/CONTCAR ./POSCAR",shell=True)
				cal2, E2_stbulk, CONTCAR_2 = second_cal("./POSCAR", log, th, ep)

				mkdir("../3_geoopt_slab")
				os.chdir("../3_geoopt_slab")
				subprocess.call('cp ../../input/* ./;  mv INCAR3 INCAR; rm INCAR1 INCAR2 INCAR4 INCAR5 POSCAR; cp ../2_constrained_bulk/CONTCAR ./POSCAR',shell=True)
				cal3, E3_slab, CONTCAR_3 = third_cal("./POSCAR", log, dlayer, nlayer, flayer, vacuum_length, slab=True)

				mkdir("../4_scf")
				os.chdir("../4_scf")
				subprocess.call('cp ../../input/* ./; mv INCAR4 INCAR; rm INCAR1 INCAR2 INCAR3 INCAR5 POSCAR; cp ../3_geoopt_slab/CONTCAR ./POSCAR',shell=True)   
				cal4, E4_scf = forth_cal(log)

				mkdir("../5_nonscf_for_wf")
				os.chdir("../5_nonscf_for_wf")
				subprocess.call('cp ../../input/* ./; mv INCAR5 INCAR; rm INCAR1 INCAR2 INCAR3 INCAR4 POSCAR; cp ../4_scf/CONTCAR ./POSCAR; cp ../4_scf/CHGCAR ../4_scf/WAVECAR ./',shell=True)
				cal5, E5_nonscf = fifth_cal(log)


			elif mode==2:
				mkdir("%s/2_constrained_bulk"%(process))
				os.chdir("%s/2_constrained_bulk"%(process))
				subprocess.call("cp ../../input/* ./; mv INCAR2 INCAR; rm INCAR1 INCAR3 INCAR4 INCAR5",shell=True)
				cal2, E2_stbulk, CONTCAR_2 = second_cal("./POSCAR", log, th, ep)

				mkdir("../3_geoopt_slab")
				os.chdir("../3_geoopt_slab")
				subprocess.call('cp ../../input/* ./;  mv INCAR3 INCAR; rm INCAR1 INCAR2 INCAR4 INCAR5 POSCAR; cp ../2_constrained_bulk/CONTCAR ./POSCAR',shell=True)
				cal3, E3_slab, CONTCAR_3 = third_cal("./POSCAR", log, dlayer, nlayer, flayer, vacuum_length, slab=True)

				mkdir("../4_scf")
				os.chdir("../4_scf")
				subprocess.call('cp ../../input/* ./; mv INCAR4 INCAR; rm INCAR1 INCAR2 INCAR3 INCAR5 POSCAR; cp ../3_geoopt_slab/CONTCAR ./POSCAR',shell=True)   
				cal4, E4_scf = forth_cal(log)

				mkdir("../5_nonscf_for_wf")
				os.chdir("../5_nonscf_for_wf")
				subprocess.call('cp ../../input/* ./; mv INCAR5 INCAR; rm INCAR1 INCAR2 INCAR3 INCAR4 POSCAR; cp ../4_scf/CONTCAR ./POSCAR; cp ../4_scf/CHGCAR ../4_scf/WAVECAR ./',shell=True)
				cal5, E5_nonscf = fifth_cal(log)

			elif mode==3:
				mkdir("%s/3_geoopt_slab"%(process))
				os.chdir("%s/3_geoopt_slab"%(process))
				subprocess.call('cp ../../input/* ./;  mv INCAR3 INCAR; rm INCAR1 INCAR2 INCAR4 INCAR5;',shell=True)
				cal3, E3_slab, CONTCAR_3 = third_cal("./POSCAR", log, dlayer, nlayer, flayer, vacuum_length, slab=True)

				mkdir("../4_scf")
				os.chdir("../4_scf")
				subprocess.call('cp ../../input/* ./; mv INCAR4 INCAR; rm INCAR1 INCAR2 INCAR3 INCAR5 POSCAR; cp ../3_geoopt_slab/CONTCAR ./POSCAR',shell=True)   
				cal4, E4_scf = forth_cal(log)

				mkdir("../5_nonscf_for_wf")
				os.chdir("../5_nonscf_for_wf")
				subprocess.call('cp ../../input/* ./; mv INCAR5 INCAR; rm INCAR1 INCAR2 INCAR3 INCAR4 POSCAR; cp ../4_scf/CONTCAR ./POSCAR; cp ../4_scf/CHGCAR ../4_scf/WAVECAR ./',shell=True)
				cal5, E5_nonscf = fifth_cal(log)

			elif mode==4:
				mkdir("%s/4_scf"%(process))
				os.chdir("%s/4_scf"%(process))
				subprocess.call('cp ../../input/* ./; mv INCAR4 INCAR; rm INCAR1 INCAR2 INCAR3 INCAR5',shell=True)       
				cal4, E4_scf = forth_cal(log)

				mkdir("../5_nonscf_for_wf")
				os.chdir("../5_nonscf_for_wf")
				subprocess.call('cp ../../input/* ./; mv INCAR5 INCAR; rm INCAR1 INCAR2 INCAR3 INCAR4 POSCAR; cp ../4_scf/CONTCAR ./POSCAR; cp ../4_scf/CHGCAR ../4_scf/WAVECAR ./',shell=True)
				cal5, E5_nonscf = fifth_cal(log)

			elif mdoe==5:
				# CHGCAR and WAVECAR must involved in 'input' directory.
				mkdir("%s/5_nonscf_for_wf"%(process))
				os.chdir("%s/5_nonscf_for_wf"%(process))
				subprocess.call('cp ../../input/* ./; mv INCAR5 INCAR; rm INCAR1 INCAR2 INCAR3 INCAR4',shell=True)
				cal5, E5_nonscf = fifth_cal(log)

			os.chdir("../../")				
			log.close()


			log=open("%s/[sss]calulation_log"%(process),'a')
			print("\n-----------Electronic properties",file=log,end="\n\n")
			print("|____1st : Surface free energy___________",file=log)
			try:
				Area=np.outer(CONTCAR_3[2][0], CONTCAR_3[2][1])
				natoms_slab=sum(CONTCAR_3[4])
				
				SFE=(E3_slab/sum(CONTCAR_3[4])-E1_bulk/sum(CONTCAR_1[4]))*natoms_slab/(2*Area)
				print(">>> Surface free energy : [ %.3f ]"%SFE, file=log)
			except:
				print(">>> Not supported yet. I will add the manual version. For now, please use calculation mode 1 or 2 to get this value", file=log)

			log.close()



	## property calculation. Later, all data files will be wrapped in .json files. (will be updated)
	# Area=np.outer(which a, which b)
	# 
	# Surface_Tension=(E3_slab-E2_stbulk)/(2*Area)
	# print("Surface tension value is [ %10.5f ], (eV/Angstrom^2)"%(Surface_Tension), file=log)
	# execute gyp.py