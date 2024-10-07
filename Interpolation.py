# -*- coding: utf-8 -*-
"""
Created the 2024.10.07
@author: Joé Petrazoller, LEM3, Metz, France : joe.petrazoller@univ-lorraine.fr
"""

#===========================================
#ONLY THE LAST TIMESTEP SHOULD APPEAR IN THE LAMMPS DUMP FILE
#Plane strain assumption : no variation in Z
#Set all atom type to 1
#The (0,0,0) position should be at the bottom left of the box
#Interpolation on a 0.1 Angstrom grid step
#===========================================


####################  Parameters to change  ####################
sigma=1 #Sigma parameter of the gaussian, in Anstrom
cutoff=3*sigma #defines the cutoff distance in Angstrom
step=0.1 #(Angstrom)
size_X=5.8849960090920456e+02 #Size of the box in the normale direction to the GB plane, in Angstrom
size_Y=5.8841678485105156e+01 #Size of the box along the GB plane, in Angstrom
filename='RESULTS.atom'     #LAMMPS dump file contraing only the last increment
################################################################

import os
import glob
import numpy as np
import pandas as pd
from math import exp



def reading_file(line):
    atom=[]
    for i in range (0,len(line)):
        tempo=[]
        tempo.append(i+1)
        atom.append(tempo)                             #"atom" is to setup a "list of list" (else we can't append after a float in the next sections)
    return atom
    
def assignin(line):
    #Creating an "atom" list containing, on each line, a list with : the atom number, the atom position x, y and z, the 6 stress components and the voronoi volume of the atom
    atom=reading_file(line)             
    for i in range(0,len(line)):
        split=line[i].split()
        number=int(split[0])
        # type_atom=split[1]
        x=float(split[2])
        y=float(split[3])
        z=float(split[4])
        sxx=float(split[5])
        syy=float(split[6])
        szz=float(split[7])
        sxy=float(split[8])
        sxz=float(split[9])
        syz=float(split[10])
        voro = float(split[11])
        atom[number-1].append(x)
        atom[number-1].append(y)  
        atom[number-1].append(z)
        atom[number-1].append(sxx*0.1/voro)  
        atom[number-1].append(syy*0.1/voro)
        atom[number-1].append(szz*0.1/voro)  
        atom[number-1].append(sxy*0.1/voro)
        atom[number-1].append(sxz*0.1/voro)  
        atom[number-1].append(syz*0.1/voro)
        atom[number-1].append(voro)          
    df1=pd.DataFrame(atom)
    max_x=max(df1[1])
    max_y=max(df1[2])
    return max_x,max_y,atom
    
def weighted_stress(x,y,max_x,atom,j,max_y,copie_up,copie_down,copie_left,copie_right,copie_topright,copie_topleft,copie_bottomright,copie_bottomleft,cutoff):
    #Interpolation core
    stress=[]
    poids=[]
    weighted_stress=[]
    for i in range(0,len(atom)):
# 9 scenario : normal box and : up, down, left, right, top right, top left, bottom right and bottom left shifts
           if y<cutoff:
               if x<cutoff: #Bottom left
                   if copie_bottomleft[i][1]-x<cutoff and copie_bottomleft[i][2]-y<cutoff:
                    distance=((copie_bottomleft[i][1]-x)**2+(copie_bottomleft[i][2]-y)**2)**0.5   #((xtarget-actualx)^2+(ytarget-actualy)^2)^0.5
                    if distance<cutoff:
                        stress.append(float(copie_bottomleft[i][4+j]))   #4+j because the 4 (5th position) of atom_window is the Sxx, 4+1 is Syy, etc
                        Wi=exp((-distance**2)/(2*(sigma**2)))         #The weight, the gaussian kernel
                        poids.append(float(Wi))
               elif x>(max_x-cutoff): #Bottom right
                   if copie_bottomright[i][1]-x<cutoff and copie_bottomright[i][2]-y<cutoff:
                    distance=((copie_bottomright[i][1]-x)**2+(copie_bottomright[i][2]-y)**2)**0.5   
                    if distance<cutoff:
                        stress.append(float(copie_bottomright[i][4+j]))   
                        Wi=exp((-distance**2)/(2*(sigma**2)))      
                        poids.append(float(Wi))
               else: #Bottom
                  if copie_down[i][1]-x<cutoff and copie_down[i][2]-y<cutoff:
                   distance=((copie_down[i][1]-x)**2+(copie_down[i][2]-y)**2)**0.5   
                   if distance<cutoff:
                       stress.append(float(copie_down[i][4+j]))   
                       Wi=exp((-distance**2)/(2*(sigma**2)))    
                       poids.append(float(Wi))
           elif y>(max_y-cutoff):
               if x<cutoff: #Top left
                   if copie_topleft[i][1]-x<cutoff and copie_topleft[i][2]-y<cutoff:
                    distance=((copie_topleft[i][1]-x)**2+(copie_topleft[i][2]-y)**2)**0.5  
                    if distance<cutoff:
                        stress.append(float(copie_topleft[i][4+j]))   
                        Wi=exp((-distance**2)/(2*(sigma**2)))         
                        poids.append(float(Wi))
               elif x>(max_x-cutoff):#Top right
                   if copie_topright[i][1]-x<cutoff and copie_topright[i][2]-y<cutoff:
                    distance=((copie_topright[i][1]-x)**2+(copie_topright[i][2]-y)**2)**0.5   
                    if distance<cutoff:
                        stress.append(float(copie_topright[i][4+j]))  
                        Wi=exp((-distance**2)/(2*(sigma**2)))         
                        poids.append(float(Wi))
               else: #Top
                   if copie_up[i][1]-x<cutoff and copie_up[i][2]-y<cutoff:
                    distance=((copie_up[i][1]-x)**2+(copie_up[i][2]-y)**2)**0.5   
                    if distance<cutoff:
                        stress.append(float(copie_up[i][4+j]))   
                        Wi=exp((-distance**2)/(2*(sigma**2)))         
                        poids.append(float(Wi))
           else:
               if x<cutoff: #Left
                   if copie_left[i][1]-x<cutoff and copie_left[i][2]-y<cutoff:
                    distance=((copie_left[i][1]-x)**2+(copie_left[i][2]-y)**2)**0.5  
                    if distance<cutoff:
                        stress.append(float(copie_left[i][4+j]))  
                        Wi=exp((-distance**2)/(2*(sigma**2)))        
                        poids.append(float(Wi))
               elif x>(max_x-cutoff): #Right
                   if copie_right[i][1]-x<cutoff and copie_right[i][2]-y<cutoff:
                    distance=((copie_right[i][1]-x)**2+(copie_right[i][2]-y)**2)**0.5   
                    if distance<cutoff:
                        stress.append(float(copie_right[i][4+j]))  
                        Wi=exp((-distance**2)/(2*(sigma**2)))         
                        poids.append(float(Wi))
               if atom[i][1]-x<cutoff and atom[i][2]-y<cutoff: #Normal configuration
                distance=((atom[i][1]-x)**2+(atom[i][2]-y)**2)**0.5  
                if distance<cutoff:
                    stress.append(float(atom[i][4+j]))   
                    Wi=exp((-distance**2)/(2*(sigma**2)))     
                    poids.append(float(Wi))
    for p in range(0,len(stress)):
        weighted_stress.append(poids[p]*stress[p])      #Weight * stress
    weighted_average=sum(weighted_stress)/sum(poids)    #Sum(Weight*stress)/Sum(Weight)
    return weighted_average, cutoff
    
def calculation(size_X,size_Y,max_x,max_y,interpolation_matrix,j,filename,atom,copie_up,copie_down,copie_left,copie_right,copie_topright,copie_topleft,copie_bottomright,copie_bottomleft,cutoff,step):
    if j==0:
        NAME='XX'
    if j==1:
        NAME='YY'
    if j==2:
        NAME='ZZ'
    if j==3:
        NAME='XY'
    if j==4:
        NAME='XZ'
    if j==5:
        NAME='YZ'
# ===== in A ======
    x_loop=size_X
    y_loop=size_Y
 # ===========   
    x_loop=x_loop*(1/step) #Ex : A 500 A long box with a 0.1 step leads to 500*10 = 5000 values
    y_loop=y_loop*(1/step)
#======================  Main loop, sweaping all the box =========================
    for x in range(0,int(10)):
        x1=x*step
        print('________'+NAME+'   ---   '+str(x) + '(it.number)')
        print('x1= '+'{:.1f}'.format(x1)+'/'+str(int(x_loop*step)))
        for y in range(0,int(10)):
            y1=y*step
            average,cutoff=weighted_stress(x1,y1,max_x,atom,j,max_y,copie_up,copie_down,copie_left,copie_right,copie_topright,copie_topleft,copie_bottomright,copie_bottomleft,cutoff)
            interpolation_matrix[y][x]=average #putting the result from "weighted_stress" in the right cell of the "interpolation_matrix" created
#=================================================================================
    # interpolation_matrix=np.flipud(interpolation_matrix)                                            # to flip the y data
    df=pd.DataFrame(interpolation_matrix)
    df.to_csv(NAME+str(cutoff)+'_sigma_'+str(sigma)+'.csv', index=False, line_terminator='\r\n')

def main():
    global atom
    global atompd
    global line
    global copie_uppd
    global sigma
    WORKDIR = glob.glob(os. getcwd())[0] #Define the Working Directory (the current directory)
    with open(filename,'r') as file:
        line=file.readlines()
        del line[0:9] #To keep only the atom positions and remove the first 9 lines header
#Define all the configurations to deal with periodic boundary conditions
    max_x,max_y,atom=assignin(line)
    max_x,max_y,copie_up=assignin(line)
    max_x,max_y,copie_down=assignin(line)
    max_x,max_y,copie_topright=assignin(line)
    max_x,max_y,copie_bottomright=assignin(line)
    max_x,max_y,copie_topleft=assignin(line)
    max_x,max_y,copie_bottomleft=assignin(line)
    max_x,max_y,copie_left=assignin(line)
    max_x,max_y,copie_right=assignin(line)

    for i in range(0,len(atom)):
        #Create "bottom" configurations      
        if atom[i][2]>(max_y-cutoff):
            copie_down[i][2]=atom[i][2]-size_Y
            copie_bottomleft[i][2]=atom[i][2]-size_Y
            copie_bottomright[i][2]=atom[i][2]-size_Y 
        if atom[i][1]>(max_x-cutoff):
            copie_bottomleft[i][1]=atom[i][1]-size_X
            copie_left[i][1]=atom[i][1]-size_X
        if atom[i][1]<(cutoff):
            copie_bottomright[i][1]=atom[i][1]+size_X     
            copie_right[i][1]=atom[i][1]+size_X
        #Create "top" configurations   
        if atom[i][2]<cutoff:
            copie_up[i][2]=atom[i][2]+size_Y
            copie_topleft[i][2]=atom[i][2]+size_Y
            copie_topright[i][2]=atom[i][2]+size_Y
        if atom[i][1]>(max_x-cutoff):
            copie_topleft[i][1]=atom[i][1]-size_X    
            copie_left[i][1]=atom[i][1]-size_X
        if atom[i][1]<(cutoff):
            copie_topright[i][1]=atom[i][1]+size_X       
            copie_right[i][1]=atom[i][1]+size_X

# ======================     To export a .csv file of the modified configurations =====================
    # copie_bottomleftpd=pd.DataFrame(copie_bottomleft)
    # copie_bottomrightpd=pd.DataFrame(copie_bottomright)
    # copie_topleftpd=pd.DataFrame(copie_topleft)
    # copie_toprightpd=pd.DataFrame(copie_topright)
    # atompd=pd.DataFrame(atom)
    # copie_uppd=pd.DataFrame(copie_up)
    # copie_downpd=pd.DataFrame(copie_down)
    # copie_leftpd=pd.DataFrame(copie_left)
    # copie_rightpd=pd.DataFrame(copie_right)   

    # atompd.to_csv('atompd.csv', index=False, line_terminator='\r\n')
    # copie_bottomleftpd.to_csv('copie_bottomleftpd.csv', index=False, line_terminator='\r\n')
    # copie_bottomrightpd.to_csv('copie_bottomrightpd.csv', index=False, line_terminator='\r\n')
    # copie_topleftpd.to_csv('copie_topleftpd.csv', index=False, line_terminator='\r\n')
    # copie_toprightpd.to_csv('copie_toprightpd.csv', index=False, line_terminator='\r\n')
    # copie_uppd.to_csv('copie_uppd.csv', index=False, line_terminator='\r\n')
    # copie_downpd.to_csv('copie_downpd.csv', index=False, line_terminator='\r\n')
    # copie_leftpd.to_csv('copie_leftpd.csv', index=False, line_terminator='\r\n')
    # copie_rightpd.to_csv('copie_rightpd.csv', index=False, line_terminator='\r\n')
#================================================================================================

 #for j in range(0,6) to treat all the 6 symmetric tensor values
    for j in range(0,6):
        #Creates the matrix that will be filled with the interpolated values, on a (Size_X * 1/step, Size_Y * 1/step grid)
        interpolation_matrix = np.zeros((int(size_Y*(1/step)),int(size_X*(1/step))),dtype=object)
        #Starting the interpolation
        calculation(size_X,size_Y,max_x,max_y,interpolation_matrix,j,filename,atom,copie_up,copie_down,copie_left,copie_right,copie_topright,copie_topleft,copie_bottomright,copie_bottomleft,cutoff,step)
        
if __name__ == '__main__':  
    main()





