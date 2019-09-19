#! /usr/bin/env python
import macrodensity as md
import math
import numpy as np
import argparse
import os

##############################################################
## -------------------Explanation------------------
## If you use this code in server (which is not local), 'plt.show()' command will not work if you don't use 'Xwindows' in Server.
## $ ssh -X [user@server] --> connect server through this command will help you.
## And you must set 'XForwarding yes' in your local terminal. You might can find in '/etc/ssh/sshd_config' file.
## ------------------------------------------------
##############################################################


pars = argparse.ArgumentParser()
pars.add_argument('-i',type=str,default='LOCPOT',help="* Name of input file. default=LOCPOT")
pars.add_argument('-o', type=str,default='planar.dat',help="* Name of output file (data file). If fermi value is existed (both manually input or extracted), it automatically shift the value.")
pars.add_argument("-fermi", help="* Set fermi level Manually", type=float, default=0)
pars.add_argument('--v', help="* If you want to see the plot, turn this on", action='store_true')
pars.add_argument('-igor',help='* If you want to get .itx file(Igor format), turn this on', action="store_true")
pars.add_argument('-igor_output',type=str,default='PlanarAverage',help="* Name of Igor output file")
pars.add_argument("--average", help="* Display Macro Average. In this case, you need lattice vector values by setting '-lv' parameter.", action='store_true')
pars.add_argument('-lv',type=float,default=2,help="* Lattice vector along the direction")
args = pars.parse_args()
input_file, Macro_average, fermi_e, output_file, visualization, igor, igor_output=args.i, args.average, args.fermi, args.o, args.v, args.igor, args.igor_output

def write_Igor2d(x,y,output, axis_name=""):
    '''x and y is raw data (np.array type).
    output is name of output file.'''
    import os
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
    itx.write("X Label left \"\Z24\F'Times New Roman'\\f02E\\f00\BPOT\M\Z24 (eV) (\\f02E\\f00-\\f02E\\f00\Bf\M)\"\n")
    itx.write("X ModifyGraph axThick=2\n")
    itx.write("X ModifyGraph lsize=2\n")
    itx.write("X ModifyGraph lblMargin(left)=5\n")
    itx.write("X Label bottom \"Distance along \\f02z\\f00 (\\{num2char(197)})\"\n")
    itx.close()


log=open('workfunction.log','w')

if os.path.exists(input_file):
    if not(os.path.exists("OUTCAR")):
        print("OUTCAR file doesn't existed")
        log.write("OUTCAR file doesn't existed\n")


    if Macro_average:
        lattice_vector=args.lv


    if args.v:
        import matplotlib.pyplot as plt
    else:
        # for non-interactive backend.
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

    log.write("input file = "+input_file)
    log.write("\nPlot Macro Average (T/F) = "+str(Macro_average))
    log.write("\noutput file = "+output_file+'\n')

    #------------------------------------------------------------------
    # Get Fermi Energy Energy from OUTCAR
    #------------------------------------------------------------------
    def extract(keyword, file='OUTCAR'):
        with open(file, 'r') as out:
            for line in out:
                if keyword in line:
                    return line.split()

    shifting=True            
    if fermi_e==0:
        try:
            # fermi_e = float(extract('efermi','vasprun.xml')[2])
            fermi_e=float(extract('E-fermi','OUTCAR')[2])
            log.write("--------------------\n")
            log.write("Fermi energy= [ "+str(fermi_e)+" ] eV.\n")
        except:
            shifting=False
            log.write("--------------------\n")
            log.write("[WARNING] Fermi energy is not extracted. Shifting is not done\n")
            print("[WARNING] fermi energy is not extracted. Shifting is not done")
    else:
        log.write("--------------------\n")
        log.write("(Manually set) Fermi energy= ["+str(fermi_e)+" ] eV.\n")
        log.write("--------------------\n")


    #------------------------------------------------------------------
    # Get the potential
    #------------------------------------------------------------------
    vasp_pot, NGX, NGY, NGZ, Lattice = md.read_vasp_density(input_file,quiet=True)
    vector_a,vector_b,vector_c,av,bv,cv = md.matrix_2_abc(Lattice)
    resolution_x = vector_a/NGX
    resolution_y = vector_b/NGY
    resolution_z = vector_c/NGZ
    grid_pot, electrons = md.density_2_grid(vasp_pot,NGX,NGY,NGZ)

    #------------------------------------------------------------------
    # POTENTIAL
    #------------------------------------------------------------------
    planar = md.planar_average(grid_pot,NGX,NGY,NGZ)
    shifted_planar=planar-fermi_e
    plt.plot(shifted_planar)


    #------------start: MACROSCOPIC AVERAGE
    if Macro_average:
        macro  = md.macroscopic_average(planar,lattice_vector,resolution_z)
        plt.plot(macro)
    #------------end: MACROSCOPIC AVERAGE

    plt.savefig('Planar.eps')
    plt.show()
    np.savetxt(output_file,shifted_planar)



    #------------------------------------------------------------------
    # Calculate Workfunction Value. (Workfunction = Vacuum Level - Fermi Level)
    #------------------------------------------------------------------
    max_data=[0,0]
    i=0
    for a in planar:
        i+=1
        if max_data[1]<a:
            max_data=[i,a]
    log.write("Vacuum Level is [ %5.10s ] eV. (Index : %s / %s)\n"%(max_data[1], max_data[0], len(planar)))
    if shifting:
        log.write("So Workfunction is [ %5.10s ] eV.\n"%(max_data[1]-fermi_e))
    else:
        log.write("Workfunction = [ None ]\n")
    log.write("--------------------\n")

    if igor==True:
        if ".itx" not in igor_output:
            igor_output+=".itx"
        index=np.arange(1,len(planar)+1,1)
        index=index*vector_c/float(len(planar))
        write_Igor2d(index,shifted_planar,igor_output,"Ep")
        log.write("Successfully finished writing itx file which name is [ %s ]"%(igor_output))


    print("--------------Workfunction Calculation is Done Please kindly look at [ workfunction.log ] file.--------------")
    #------------------------------------------------------------------
    # Successfully Ended
    #------------------------------------------------------------------
else:
    log.write(input_file+" doesn't existed... Just Ended")
    log.close()

#------------------------------------------------------------------
# Not ended well
#------------------------------------------------------------------