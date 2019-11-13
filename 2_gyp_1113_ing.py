#! /usr/bin/env python
import numpy as np
import argparse
import os
import sys

##########################################################################
## -----------------------------About code----------------------------- ##
## Giyeok Lee, Ethan                                                    ##
## sss_package of MTG (https://github.com/materials-theory/sss_package) ##
## Version : 2019.11.05 (By Giyeok Lee)                                 ##
## Revised : 2019.11.05 (By Giyeok Lee)                                 ##
## -------------------------------------------------------------------- ##
##########################################################################

###################################################################################################################################################
## --------------------------------------------------------------------------------------------------------------------------------------------- ##
## If you use this code in server (which is not local), 'plt.show()' command will not work if you don't use 'Xwindows' in Server.                ##
## $ ssh -X [user@server] --> connect server through this command will help you.                                                                 ##
## And you must set 'XForwarding yes' in your local terminal. You might can find in '/etc/ssh/sshd_config' file.                                 ##
## If troubles using $DISPLAY keep happens, tag out all parts of Visualization (2), which is matplotlib part. It will solve that problem.        ##
## --------------------------------------------------------------------------------------------------------------------------------------------- ##
###################################################################################################################################################


pars = argparse.ArgumentParser()
pars.add_argument('-i',type=str,default='LOCPOT',help="* Name of input file. default=LOCPOT")
pars.add_argument('-o', type=str,default='output.dat',help="* Name of output file (data file). If fermi value is existed (both manually input or extracted), it automatically shift the value.")
pars.add_argument("-fermi", help="* Set fermi level Manually", type=float, default=0)
pars.add_argument('-v', help="* If you want to see the plot, turn this on", action='store_true')
pars.add_argument('-igor',help='* If you want to get .itx file(Igor format), turn this on', action="store_true")
pars.add_argument('-igor_o',type=str,default='Potential',help="* Name of Igor output file")
pars.add_argument('-d',type=str,default='z',help="Wanted direction. You can set a/b/c or x/y/z or 1/2/3 for the direction.")
pars.add_argument('--vac_Ediff', default=1E-3, type=float, help="Vacuum region convergence Ediff")
pars.add_argument('--vac_width', default=3.0, type=float, help="Vacuum region convergence width. Unit is angstrom")
args = pars.parse_args()
input_file, fermi_e, output_file, visualization, igor, igor_output,direction=args.i, args.fermi, args.o, args.v, args.igor, args.igor_o, args.d
vac_con_Ediff, vac_con_width = args.vac_Ediff, args.vac_width

def mkdir(dirname):
    if not os.path.exists(os.path.dirname(dirname+"/")):
        os.makedirs(os.path.dirname(dirname+"/"))


def read_CHGCAR(ipf="LOCPOT",direction="Z"):
    ip=open(ipf,'r')
    print ('* Name of System : ' + ip.readline())
    ip.readline()
    cell_vec=[]
    cell_vec.append(ip.readline().split())
    cell_vec.append(ip.readline().split())
    cell_vec.append(ip.readline().split())
    cell_vec=np.array(cell_vec,dtype="d")
    species=ip.readline().split()
    natoms=list(map(lambda x:int(x),ip.readline().split()))
    tot_natoms=sum(natoms)
    coord_read=ip.readline()
    if coord_read.startswith("D") or coord_read.startswith("d"):
        coord_type="Direct"
    elif coord_read.startswith("F") or coord_read.startswith("f"):
        coord_type="Direct"
    elif coord_read.startswith("C") or coord_read.startswith("c"):
        coord_type="Cartesian"
    else:
        print("------- Coordination Type Recognition Failed. Calculation ended. -------")
        sys.exit(1)

    coord=[]
    for i in range(tot_natoms):
        coord.append(list(map(lambda x:float(x),ip.readline().split())))
    ip.readline()
    grids = list(map(lambda x:int(x),ip.readline().split()))
    nx, ny, nz = grids
    print("* Matrix : [ %d ] x [ %d ] x [ %d ]"%(nx, ny, nz))
    lines = ip.readlines()

#------------------------------------------------------------
# From CHG.py (WS)
    pot = []
    for x in lines:
        if "a" in x:
            break
        else:
            tmp = x.split()
            for y in tmp:
                pot.append(y)
    pot = np.reshape(np.array(pot, dtype="d"), (grids[::-1]))
    # if you want to undo reshape, use command, pot.flatten()
#------------------------------------------------------------
    if direction in "Z":
    # since chglist is nz ny nx sequence, so if we want nz direction, mean axis will be axis=(1,2)
    # if axis is in tuple type, they calculate average in both direction.
        potavg=np.mean(pot, axis=(1,2))
    if direction in "Y":
        potavg=np.mean(pot, axis=(2,0))
    if direction in "X":
        potavg=np.mean(pot, axis=(0,1))
    ip.close()
    return [cell_vec, species, natoms, coord_type, coord], grids, pot, potavg



def write_Igor2d(x,y,lLabel, bLabel, output, axis_name=""):
    '''x and y is raw data (np.array type).
    output is name of output file.'''
    if os.path.exists(output):
        print("%s already exists. File substituted to new file."%(igor_output))
    itx=open(output,'w')
    itx.write("IGOR\n")
    itx.write("WAVES/D x_%s %s\n"%(axis_name,axis_name))
    itx.write("BEGIN\n")
    for (x1,y1) in zip(x,y):
        itx.write(str(x1)+" "+str(y1)+"\n")
    itx.write("END\n")
    itx.write("X Display %s vs x_%s as \"planar_average_%s\" \n"%(axis_name,axis_name,axis_name))
    itx.write("X DefaultFont/U \"Times New Roman\"\n")
    itx.write("X ModifyGraph marker=19\n")
    itx.write("X ModifyGraph tick=2\n")
    itx.write("X ModifyGraph mirror=1\n")
    itx.write("X ModifyGraph fSize=28\n")
    itx.write("X ModifyGraph lblMargin(left)=15,lblMargin(bottom)=10\n")
    itx.write("X ModifyGraph standoff=0\n")
    itx.write("X ModifyGraph axThick=1.5\n")
    itx.write("X ModifyGraph axisOnTop=1\n")
    itx.write("X ModifyGraph width=453.543,height=340.157\n")
    itx.write("X ModifyGraph zero(left)=8\n")
    itx.write("X ModifyGraph zero(bottom)=1\n")
    itx.write("X "+lLabel)
    itx.write("X ModifyGraph axThick=2\n")
    itx.write("X ModifyGraph lsize=2\n")
    itx.write("X ModifyGraph lblMargin(left)=5\n")
    itx.write("X "+bLabel)
    itx.close()

def extract(keyword, file='OUTCAR'):
    '''Get Fermi Energy Energy from OUTCAR'''
    with open(file, 'r') as out:
        for line in out:
            if keyword in line:
                return line.split()

def length(a):
    return (a[0]**2+a[1]**2+a[2]**2)**(0.5)

def dir2car(coord,cell_vec):
    '''Direct to Cartesian (function from gycoco.py)
    [input] : coord, cell_vec
    |-> coord : coordination. @np.array(dtype="d")
    |-> cell_vec : cell vector. @np.array(dtype='d')
    [output] : new coordination @np.array(dtype='d')
    '''
    return np.matmul(coord,cell_vec)


#------------------------------------------------------------------
# Starts Script.
#-----------------------------------------------------------------
#------------------------------------------------------------------
# input Direction testing
#------------------------------------------------------------------

if direction in "zZcC3":
    direction="Z"
    ndir=3
elif direction in "yYbB2":
    direction="Y"
    ndir=2
elif direction in "xXaA1":
    direction="X"
    ndir=1
else:
    raise IOError("------- Wrong input of direction. Please set direction again -------")
    sys.exit(1)

#------------------------------------------------------------------
# All files are generated in this directory.
#------------------------------------------------------------------
directory="wf_cal"
mkdir(directory)
log=open(directory+'/workfunction.log','w')

#------------------------------------------------------------------
# Extract fermi energy if 'OUTCAR' file exists in same directory
#------------------------------------------------------------------
if os.path.exists(input_file):
    print("--------------        {} exists. Calculation starts.        --------------".format(input_file))
    log.write(input_file+" exists. Calculation starts\n")
    log.write("--------------------\n")
    log.write("input file = "+input_file+'\n')
    log.write("output file = "+output_file+'\n')

    #------------------------------------------------------------------
    # Get the potential
    #-----------------------------------------------------------------
    values=read_CHGCAR(input_file, direction)
    potavg, grids, cell_vec = values[3], np.array(values[1],dtype='d'), np.array(values[0][0],dtype='d')
    cell_vec_length = np.array((length(values[0][0][0]), length(values[0][0][1]), length(values[0][0][2])),dtype='d')
    resolution=cell_vec_length/grids

    #------------------------------------------------------------------
    # Get E_vacuum from average potential along the direction.
    #------------------------------------------------------------------
    potavg_temp=np.append(np.array([0],dtype='d'),potavg[:-1])
    diff_potavg=abs(potavg-potavg_temp)


    log.write("--------------------\n")
    print(">>> Vacuum Level Searching",file=log)
    print("|---> E_diff for vacuum convergence is %f"%(vac_con_Ediff),file=log)


    temp, x_temp =[], []
    for i in range(len(diff_potavg)):
        if diff_potavg[i]<(vac_con_Ediff):
            temp.append(potavg[i])
            x_temp.append(i)
    x_temp=np.array(x_temp,dtype='d')
    coord=np.array(values[0][4],dtype='d')

    if len(x_temp)==0:
        print("|---> No Flat region was found. Maybe Dipole correction or Increasing vacuum level is needed.", file=log)
    else:
        diff_x_temp=x_temp-np.append(np.array([0],dtype='d'),x_temp[:-1])
        n=0
        for i in diff_x_temp[1:]:
            if i!=1:
                n+=1
    if n==0:
        E_vac=np.mean(temp)
        if len(x_temp)>=vac_con_width/resolution[ndir-1]:
            print("|---> Sufficient vacuum region about %5.3f angstrom width was found."%(len(x_temp)*resolution[ndir-1]),file=log)
        else: 
            print("|---> Insufficient vacuum region about %5.3f angstrom. Temporarily, this narrow flat region was used as vacuum level. Please check about it."%(len(x_temp)*resolution[ndir-1]),file=log)
    elif n==1:
        for i in range(len(diff_x_temp[1:])):
            if diff_x_temp[i]!=1:
                discrete=i
        region_a=np.mean(x_temp[:discrete])
        region_b=np.mean(x_temp[discrete:])
        if values[0][3]=="Direct":
            x_model=np.mean(dir2car(coord,cell_vec))/resolution[ndir-1]
        else:
            x_model=np.mean(coord[:,2])/resolution[ndir-1]
        if abs(x_model-region_a) <= abs(x_model-region_b):
            region=np.array(x_temp[:discrete],dtype='d')
            E_region=np.array(temp[:discrete],dtype='d')
        else:
            region=np.array(x_temp[discrete:],dtype='d')
            E_region=np.array(temp[discrete:],dtype='d')
        E_vac=np.mean(E_region)
        print("|---> 2 Flat regions were found. Maybe dipole correction was performed. Please check about it.", file=log)
        print("|---> If you turn on visualization by -v tags, you can see blue lines. That is the selected region.", file=log)
    else:
        E_vac=max(temp)
        print("|---> Too many flat regions were found. Please check the convergence of Vacuum. For now, max energy is used.", file=log)


    #------------------------------------------------------------------
    # Shifting potential values with respect to E_fermi (default E_fermi = 0.0)
    #------------------------------------------------------------------
    shifting=True            
    if fermi_e+1==1:
        try:
            # fermi_e = float(extract('efermi','vasprun.xml')[2])
            fermi_e=float(extract('E-fermi','OUTCAR')[2])
            log.write("--------------------\n")
            log.write("Fermi energy= [ "+str(fermi_e)+" ] eV.\n")
            lLabel="Label left \"\Z24\F'Times New Roman'\\f02E\\f00\BPOT\M\Z24 (eV) (\\f02E\\f00-\\f02E\\f00\Bf\M)\"\n"
        except:
            shifting=False
            log.write("--------------------\n")
            log.write("OUTCAR file doesn't existed\n")
            log.write("[WARNING] Fermi energy is not extracted. Shifting is not done\n")
            print("OUTCAR file doesn't existed")
            print("[WARNING] fermi energy is not extracted. Shifting is not done")
            lLabel="Label left \"\Z24\F'Times New Roman'\\f02E\\f00\BPOT\M\Z24 (eV)\"\n"
    else:
        log.write("--------------------\n")
        print("OUTCAR file doesn't existed")
        print("(Manual setting) Fermi Energy is successfully set.")
        log.write("OUTCAR file doesn't existed\n")
        log.write("(Manually set) Fermi energy= ["+str(fermi_e)+" ] eV.\n")
        log.write("--------------------\n")
        lLabel="Label left \"\Z24\F'Times New Roman'\\f02E\\f00\BPOT\M\Z24 (eV) (\\f02E\\f00-\\f02E\\f00\Bf\M)\"\n"
    shifted_potavg=potavg-fermi_e

    #------------------------------------------------------------------
    # Calculate Workfunction Value. (Workfunction = Vacuum Level - Fermi Level)
    #------------------------------------------------------------------

    log.write("Vacuum Level is [ %5.10s ] eV.\n"%(E_vac))
    if shifting:
        log.write("So Workfunction is [ %5.10s ] eV.\n"%(E_vac-fermi_e))
    else:
        log.write("Workfunction = [ None ]\n")
    log.write("--------------------\n")



    #------------------------------------------------------------------
    # Visualization -- Total
    #------------------------------------------------------------------
    index=np.arange(1,len(potavg)+1,1)
    if direction == "Z":
        real_index=index*cell_vec_length[2]/float(len(potavg))
    elif direction == "Y" :
        real_index=index*cell_vec_length[1]/float(len(potavg))
    else:
        # direction == "X"
        real_index=index*cell_vec_length[0]/float(len(potavg))

    #------------------------------------------------------------------
    # Successfully Ended.
    #------------------------------------------------------------------    
    print("--------------         Workfunction Calculation is Done.        --------------")
    print("-------------- Please kindly look at [ workfunction.log ] file. --------------")
    outdat=open(directory+"/"+output_file,'w')
    for (x1,y1) in zip(index,shifted_potavg):
        print((str(x1)+" "+str(y1)),file=outdat)
    outdat.close()
       
    
    #------------------------------------------------------------------
    # Visualization (1) Igor
    #------------------------------------------------------------------
    if igor==True:
        if ".itx" not in igor_output:
            igor_output+=".itx"
        #lLabel is written in upper parts. Since it must be differnent case by case.
        bLabel="Label bottom \"Distance along \\f02%s\\f00 (\\{num2char(197)})\"\n"%(direction.lower())
        write_Igor2d(real_index,shifted_potavg,lLabel,bLabel,directory+"/"+igor_output,"Ep")
        log.write("Successfully finished writing itx file which name is [ %s ]"%(igor_output))
        
    #------------------------------------------------------------------
    # Visualization (2) matplotlib.pyplot
    #------------------------------------------------------------------
    import matplotlib.pyplot as plt
    X, Y = real_index, shifted_potavg
    plt.plot(X,Y,'r')
    if n!=0:
        # If there are dipole corrections
        selected_X, selected_Y = region, E_region-fermi_e
        plt.plot(selected_X*resolution[ndir-1],selected_Y,'bo')
    plt.savefig(directory+"/"+'Potential.eps')
    if not(visualization):
        import matplotlib
        matplotlib.use("Agg")
    else:
        plt.show()
    
    log.close()

#------------------------------------------------------------------
# Not ended well
#------------------------------------------------------------------
else:
    log.write(input_file+" doesn't existed... Just Ended")
    log.close()