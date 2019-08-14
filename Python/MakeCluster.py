#!/usr/bin/env python

##################################################################################
#                               Author : Mehdi Zare                              #
#                               Date   : 8/4/2019                                #
#                               University of South Carolina                     #
#       Purpose: Make 51-atom cluster out of 111 or 0001 surfaces	         #
##################################################################################   
import numpy as np
import pandas as pd
from pandas import DataFrame, Series


#  Function for reading vasp input coordinate file (POSCAR/CONTCAR)
#  It needs one argument which needs to be in vasp format (like POSCAR/CONTCAR)
def vasp_read_poscar(filename):
    try:
        fhand=open(filename)
    except:
        print('file ', filename,' cannot be opened')
        return
    contcar=fhand.readlines()
    
    #read cell dimensions
    cell=[]
    cell.append(float(contcar[2].split()[0]))
    cell.append(float(contcar[3].split()[1]))
    cell.append(float(contcar[4].split()[2]))
    
    #read type of coordinates: cartesian or fractional
    mode=contcar[8].strip()
    
    #read information about number of atoms and types of them
    num_atom=len(contcar[5].split())
    atoms={}
    for index in range(0,num_atom):
        atoms[contcar[5].split()[index]]=int(contcar[6].split()[index])
    
    coordinates={}
    initial=9
    for key, values in atoms.items():
        position=[]
        for i in range(initial,values+initial):
            position.append([float(k) for k in contcar[i].split()[0:3] ])
            
        #print(len(position))
        initial=values+initial
        coordinates[key]=position

    return mode, cell, atoms, coordinates, contcar
    fhand.close()


# # Function for conveting POSCAR coordinates
# 
# ### This function needs two input: filename which needs to be a POSCAR format file and the other one is mode that could be F2C (Fractional to Cartesian) or  C2F (Cartesian to Fractional) 
# ### Example for calling: vasp_conver_poscar_coord(filename='CONTCAR',mode='F2C')
def vasp_conver_poscar_coord(filename,mode):
    readfile=vasp_read_poscar(filename)
    if mode == 'F2C':
        if readfile[0].lower() != 'direct':
            print('Your input file shoud be in Fractional coordinate')
            return
        else:
            cell=readfile[1]
            coordinates=readfile[3]
            for atoms, coord in coordinates.items():
                for i in range(0,len(coord)):
                    for j in range(0,3):
                        if coordinates[atoms][i][j] < 0: # for negative coordinates
                            coordinates[atoms][i][j]=(coordinates[atoms][i][j]+1.0) * cell[j]
                        else:
                            coordinates[atoms][i][j]=coordinates[atoms][i][j] * cell[j]
            fwrite=open(filename+'-'+mode, 'w')
            #write header up to before coordinates
            for i in range(0,8):
                fwrite.write(readfile[4][i])
            fwrite.write('Cartesian\n')
            #write Coordinates
            indexFT=9  #the starting index for fixed or relaxed atom
            for atoms, coord in coordinates.items():
                for i in range(0,len(coord)):
                    a=['{0:.15f}'.format(x) for x in coordinates[atoms][i]]
                    F_or_T='   '.join(readfile[4][indexFT].split()[3:6])
                    indexFT += 1
                    fwrite.write('       '+'     '.join(a)+'   '+F_or_T+'\n')
    
            fwrite.close() 
    elif mode == 'C2F':
        if readfile[0].lower() != 'cartesian':
            print('Your input file shoud be in Cartesian coordinate')
            return
        else:
            cell=readfile[1]
            coordinates=readfile[3]
            for atoms, coord in coordinates.items():
                for i in range(0,len(coord)):
                    coordinates[atoms][i][0]=coordinates[atoms][i][0] / cell[0]
                    coordinates[atoms][i][1]=coordinates[atoms][i][1] / cell[1]
                    coordinates[atoms][i][2]=coordinates[atoms][i][2] / cell[2]
            fwrite=open(filename+'-'+mode, 'w')
            #write header up to before coordinates
            for i in range(0,8):
                fwrite.write(readfile[4][i])
            fwrite.write('Direct\n')
            #write Coordinates
            indexFT=9  #the starting index for fixed or relaxed atom
            for atoms, coord in coordinates.items():
                for i in range(0,len(coord)):
                    a=['{0:.15f}'.format(x) for x in coordinates[atoms][i]]
                    F_or_T='   '.join(readfile[4][indexFT].split()[3:6])
                    indexFT += 1
                    fwrite.write('       '+'     '.join(a)+'   '+F_or_T+'\n')
    
            fwrite.close()     
    else:
        print('Please insert F2C or C2F')
        return
    
    return coordinates
    

# # Function for finding the information about your metal slab such as:
# ## number of layers, atoms in each layer, distance between each layer, and distance between atoms in X or Y direction
# 
# # Warning: This function is not trustworthy
def vasp_slab_info(file,metal):
    readfile=vasp_read_poscar(file)
    mode=readfile[0]  ; cell=readfile[1];
    atoms=readfile[2] ;
    
    if mode.lower() == 'direct':
        coordinates=vasp_conver_poscar_coord(filename=file,mode='F2C')
    else:
        coordinates=readfile[3]
    
    #find the maximum Z of metal and difine limits of each layer
    max_z_metal=0
    min_z_metal=100000
    for atoms, coord in coordinates.items():
        for i in range(0,len(coord)):
            Z_coord=coordinates[atoms][i][2]
            if atoms == metal and  Z_coord > max_z_metal:
                max_z_metal= Z_coord
            elif atoms == metal and Z_coord < min_z_metal:
                min_z_metal= Z_coord
                
    #I condider the movement of each layer is not more than 1 angestrom
    Top_layer=[];
    Bot_layer=[];
    for atoms, coord in coordinates.items():
        for i in range(0,len(coord)):
            Z_coord=coordinates[atoms][i][2]
            if atoms == metal:
                if   Z_coord > max_z_metal - 1.0:
                    Top_layer.append(coordinates[atoms][i])
                elif Z_coord < min_z_metal + 1.0:
                    Bot_layer.append(coordinates[atoms][i])
    #FLAG FOR ERROR               
    if len(Top_layer) != len(Bot_layer):
        print('ERROR ! THIS CODE CAN NOT HANDLE YOUR FILE, number of atoms in Top layer and bottom layer are different')
        return
      
    #atom in each layer
    atom_each_layer=len(Top_layer)
    
    diff=[0,0,0] # half of the differences in X,Y,Z direction between atoms
    #find number of layer
    atoms=readfile[2]
    layers= int((atoms[metal]-(len(Top_layer)*2))/len(Top_layer))  + 2  #2 is for TOP and BOT layers
    inter_layer=(max_z_metal/(layers-1)/2)
    limits_Z=[0.0+(inter_layer*l) for l in range(1,layers*2+1,2)]
    diff[2]=inter_layer
    
    #find Distance between atoms in X and Y dirction
    # I consider the movement in Y direction is not more than 1 Ang.
    Y_max_TOP=(max(Top_layer, key=lambda x: x[1])[1])
    Y_min_TOP=(min(Top_layer, key=lambda x: x[1])[1])
    X_max_TOP=(max(Top_layer, key=lambda x: x[0])[0])
    X_min_TOP=(min(Top_layer, key=lambda x: x[0])[0])
    
    Y_first=[];
    Y_last=[]
    for i in Top_layer:
        if i[1] < Y_min_TOP + 1.0:
            Y_first.append(i)
        elif i[1] > Y_max_TOP - 1.0:
            Y_last.append(i)
            
    #FLAG FOR ERROR
    if len(Y_first) != len(Y_last):
        print('ERROR: each layer has differnt number of atoms in every row, I cannot handle it')
        return

    
    rows_Y=int((len(Top_layer) - (len(Y_first)*2)) / len(Y_first))  + 2  # 2 for first and last rows
    rows_X=int(len(Top_layer) / rows_Y)
    
    #sort the first row based on X coordinates
    Y_first.sort(key=lambda x : x[0])
    
    #find limits of Y
    inter_layer=((Y_max_TOP-Y_min_TOP)/(rows_Y-1)/2)
    limits_Y=[0.0+(inter_layer*l) for l in range(1, rows_Y*2 +1 , 2)]
    diff[1]=inter_layer
    
    #find limits of X
    inter_layer=((X_max_TOP-X_min_TOP)/(rows_X-1)/2)
    limits_X=[0.0+(inter_layer*l) for l in range(1, rows_X*2 +1 , 2)]
    diff[0]=inter_layer
    
    
    #information of each layer
    each={};
    each['num']=len(Top_layer)
    each['rowX']=rows_X
    each['rowY']=rows_Y
    
    #limits of X, Y, Z directions
    limits={'X':limits_X,'Y':limits_Y,'Z':limits_Z}
    
    return layers, each, diff,  limits
    


# # Function for Creating xyz supercell out of CONTCAR/POSCAR
# ### It needs 6 arguments, input file name which needs to be in POSCAR format, X , Y, Z which are the reptition you want in each direction, output file name you want, and the mode which shows the format of your output file: I already only have xyz format mode, I might add vasp supercell format to it as well
# 
# ###### Example: vasp_build_supercell(file_IN='CONTCAR-F2C',mode='xyz',X=3,Y=3,Z=3,file_OUT='expanded.xyz')
def vasp_build_supercell(file_IN,mode,X,Y,Z,file_OUT):
    readfile=vasp_read_poscar(file_IN)
    cell=readfile[1]
    X=int(X)
    Y=int(Y)
    Z=int(Z)
    
    # Extract Cartesian coordinates 
    if readfile[0].lower() == 'direct':
        coordinates=vasp_conver_poscar_coord(filename=file_IN,mode='F2C')
    elif readfile[0].lower() == 'cartesian':
        coordinates=readfile[3]
        
    supercell={}
    sup_atoms=0 # number of atoms in supercell for header of xyz file
    for atoms, coord in coordinates.items():
        arr=np.array(coord)
            
        #expand in X direction
        final=arr
        for xdirect in range(1,X):
            new_arr=[[item[0]+(cell[0]*xdirect), item[1], item[2]] for item in arr]
            final=np.append(final,new_arr,axis=0)
            
        #expand in Y direction
        arr=final
        for ydirect in range(1,Y):
            new_arr=[[ item[0], item[1]+(cell[1]*ydirect), item[2]] for item in arr]
            final=np.append(final,new_arr,axis=0)
            
        #expand in Z direction
        arr=final
        for zdirect in range(1,Z):
            new_arr=[[ item[0], item[1], item[2]+(cell[2]*zdirect)] for item in arr]
            final=np.append(final,new_arr,axis=0)
            
        sup_atoms += np.shape(final)[0]
        supercell[atoms]=final
        
    #write xyz file
    #write header
    fwrite=open(file_OUT,'w')
    fwrite.write(str(sup_atoms) + '\n')
    fwrite.write('Surface' + '\n')
    #write coordinates
    for key, values in supercell.items():
        for i in range(0,len(values)):
            a=['{0:.15f}'.format(x) for x in supercell[key][i]]
            fwrite.write(key+'        '+'     '.join(a)+'\n')
    fwrite.close()
    
    return supercell
        


# # Create Cluster out of a CONTCAR/POSCAR file
# ### This gets the name of your input file, which is vasp poscar formatted, and give you a cluster in .xyz format. metal keyword is your metal name, mode can be eSMS or iSMS, and finally dX,dY,dZ are the distances between your atoms in X,Y and Z direction that you can get them easyly from your input file
# 
# # warning: this function can only create a specific shape cluster tha consists of 51 metal atom for 111 and 0001 cleaved surfaces
# 
# ### Example in eSMS: Make_Cluster_51(IN='CONTCAR', OUT='cluster.xyz',metal='Ni',X=4,Y=5,mode='esms',dX=2.5,dY=2.16,dZ=1.996), It is worth mentioning that in this mode you need to define the 5 middle atoms of cluster that you want. To get them, run this program with and it will write for you file removed_2_layers.xyz. From that file you can get those numbers
# 
# ### example in iSMS Mode: Make_Cluster_51(IN='CONTCAR', OUT='cluster.xyz',metal='Ni',X=7,Y=7,mode='isms',dX=2.5,dY=2.16,dZ=1.996)
def Make_Cluster_51(IN,OUT,X,Y,metal,mode,dX,dY,dZ):
    
    
    def dist(atom1,atom2):
        distance=np.sqrt((atom2[0]-atom1[0])**2+(atom2[1]-atom1[1])**2)
        return distance

    if mode.lower() == 'isms':
        X=7; Y=7
    
    supercell=vasp_build_supercell(file_IN=IN,mode='xyz',X=X,Y=Y,Z=1,file_OUT='expanded.xyz')    
    X=int(X)
    Y=int(Y)
    #Removing the first two layers
    
    #find the maximum Z of metal and difine limits of each layer
    max_z_metal=0
    for atoms, coord in supercell.items():
        for i in range(0,len(coord)):
            if atoms == metal and supercell[atoms][i][2] > max_z_metal:
                max_z_metal= supercell[atoms][i][2]
    inter_layer=(max_z_metal/3/2)
    limits=[0.0+(inter_layer*l) for l in range(1,9,2)]
    #print('limits=',limits)
    
    #remove the first two layer of metals based on limits and write the new surface to file removed_2_layers
    removed_2_layers=supercell
    metal_list=[]
    for atoms, coord in supercell.items():
        for i in range(0,len(coord)):
            if atoms == metal and supercell[atoms][i][2] > limits[1]:
                metal_list.append(supercell[atoms][i])
    arr_metal=np.array(metal_list)
    removed_2_layers[metal]=arr_metal
    
    ####################### write removed_2_layers #################
    tot_atoms=0
    for atoms, coord in removed_2_layers.items():
        tot_atoms += np.shape(coord)[0]
        
    
    #write xyz file
    #write header
    fwrite=open('removed_2_layers.xyz','w')
    fwrite.write(str(tot_atoms) + '\n')
    fwrite.write('Surface' + '\n')
    #write coordinates
    for key, values in removed_2_layers.items():
        for i in range(0,len(values)):
            a=['{0:.15f}'.format(x) for x in removed_2_layers[key][i]]
            fwrite.write(key+'        '+'      '.join(a)+'\n')
    fwrite.close()
    ############################ END OF write   ####################
    
    
    if mode.lower() == 'esms':
        one,two,three,four,five=[int(x) for x in input("You need to insert five numbers four me:\n"    
                            "the first two values are your middle cluster atoms in"
                            " X direction\n"
                            " the second two valuse are your middle cluster atoms in"
                            " Y direction\n, and the last value is the middle cluster"
                            " atom number from bottom layer exactl under those four"
                            " atoms you just entered,\n write them down with space").split()]
        # locate cell number
        metal_each_cell=np.shape(removed_2_layers[metal])[0]/(X*Y)
        [Q, R]=[int(one//metal_each_cell),int(one%metal_each_cell)]
        cell_num=Q if R == 0 else Q+1 
        #print(cell_num)
        
    elif mode.lower() == 'isms':
        cell_num=int(25)

    #remove all adsorbate except of cell_num
    for atoms, coord in removed_2_layers.items():
        if atoms != metal:
            atom_each_cell=int(np.shape(coord)[0]/(X*Y))
            #print(atom_each_cell)
            #Extract adsrobate of cell_num
            first_index=atom_each_cell*(cell_num - 1)
            atom_ads_list=[]
            for i in range(first_index,first_index+atom_each_cell):
                         atom_ads_list.append(removed_2_layers[atoms][i])
            arr_atom=np.array(atom_ads_list)
            removed_2_layers[atoms]=arr_atom
    
    if mode.lower() == 'isms':
        #find the box dimension of adsorbate
        Adsorbate=[]
        for atoms, coord in removed_2_layers.items():
            max_X=0;
            max_Y=0;
            if atoms != metal:
                for i in range(0,len(coord)):
                    Adsorbate.append(removed_2_layers[atoms][i])
                    #X_coord=removed_2_layers[atoms][i][0]
                    #Y_coord=removed_2_layers[atoms][i][1]
                    #if  X_coord > max_X:
                    #    max_X=X_coord
                    #elif X_coord < min_X:
                     #   min_X=X_coord
                    #if Y_coord > max_Y:
                     #   max_Y=Y_coord
                    #elif Y_coord < min_Y:
                      #  min_Y=Y_coord
        #print('AdsCoord',Adsorbate)
        max_X=max(i[0] for i in Adsorbate)
        min_X=min(i[0] for i in Adsorbate)
        max_Y=max(i[1] for i in Adsorbate)
        min_Y=min(i[1] for i in Adsorbate)
                        
        #print('adsboxX=',max_X,min_X)
        #print('adsbixY=',max_Y,min_Y)
        
        #middle of the box of adsorbate
        X_mid=min_X + (max_X-min_X)/2; Y_mid=min_Y+ (max_Y-min_Y)/2
        center=np.array([X_mid,Y_mid])
        
        #print('Xmid=',X_mid,'Ymid=',Y_mid)
        
        #half distance of atoms in X,Y,Z direction
        X_diff=dX/2
        Y_diff=dY/2
        Z_diff=dZ/2
        
        #print('Xdiff=',X_diff,'Ydiff=',Y_diff,'Zdiff=',Z_diff)
        
        #find metal atoms located in adsorbate box and a little more( 2 Angstrom shift)
        near_ads=[]
        for i in removed_2_layers[metal]:
            X_coord=i[0]
            Y_coord=i[1]
            Z_coord=i[2]
            if (Z_coord > limits[2]) and ( min_X - 2.0  < X_coord < max_X + 2.0) and (min_Y-2.0 < Y_coord < max_Y+2.0):
                near_ads.append(i)
        #print('atoms of top layer close to adsorbate',near_ads)
        
        
        #consider each of above atoms in near_ads as one and get the other three around it
        min_distance=100000.0
        for k in near_ads:
            X_one=k[0];Y_one=k[1];Z_one=k[2]
            oneCoord=k
            distance = dist(k,center)
            for i in removed_2_layers[metal]:
                X_coord=i[0];Y_coord=i[1];Z_coord=i[2]
                if Z_coord > Z_one - Z_diff: # top layer
                    if  X_one + X_diff < X_coord < X_one + 3*X_diff and Y_one - Y_diff < Y_coord < Y_one + Y_diff: #right atom
                        twoCoord=i
                        distance += dist(i,center)
                    if X_one < X_coord < X_one + 2*X_diff and Y_one - 3*Y_diff < Y_coord < Y_one - Y_diff:        #down right atom
                        threeCoord=i
                        distance += dist(i,center)
                    if X_one < X_coord < X_one + 2*X_diff and Y_one + Y_diff < Y_coord < Y_one + 3*Y_diff:        #top right atom
                        fourCoord=i
                        distance += dist(i,center)
            if distance < min_distance:
                min_distance=distance
                min_four=[oneCoord,twoCoord,threeCoord,fourCoord]
                        
        #print('\n mindistance=',min_distance)           
        #print('\n four atoms are', min_four)
        
        # find atom five
        X_one=min_four[0][0]; Y_one=min_four[0][1]; Z_one=min_four[0][2]
        for i in removed_2_layers[metal]:
            X_coord=i[0];Y_coord=i[1];Z_coord=i[2]
            if Z_coord < Z_one - Z_diff: # bottom layer
                if  X_one < X_coord < X_one + 2*X_diff and Y_one - Y_diff < Y_coord < Y_one + Y_diff:
                    min_four.append(i)
        
        
        five_atoms=min_four
        
        #print(five_atoms)
        
        [X1,Y1,Z1]=five_atoms[0]
        [X2,Y2,Z2]=five_atoms[1]
        [X3,Y3,Z3]=five_atoms[2]
        [X4,Y4,Z4]=five_atoms[3]
        [X1BOT,Y1BOT,Z1BOT]=five_atoms[4]
    
            
        
     #temperary for isms   
    #one=int(798);two=int(791);three=int(563);four=int(795);five=int(779)
    elif mode.lower() == 'esms':
        [X1,Y1,Z1]=removed_2_layers[metal][one-1]
        [X2,Y2,Z2]=removed_2_layers[metal][two-1]
        [X3,Y3,Z3]=removed_2_layers[metal][three-1]
        [X4,Y4,Z4]=removed_2_layers[metal][four-1]
        [X1BOT,Y1BOT,Z1BOT]=removed_2_layers[metal][five-1]
        
        
    Xdiff=(X2-X1)/2; Ydiff=(Y4-Y1)/2; Xcenter=X1+Xdiff
    
    #bonds for the middle row top layer
    XpRP0=X2+(5*Xdiff);XmRP0=X1-(5*Xdiff)
    Yp1=Y1+Ydiff; Ym1=Y1-Ydiff
    
    #Row Plus1
    XpRP1=X2+(4*Xdiff);XmRP1=X1-(4*Xdiff)
    Yp2=Y1+(3*Ydiff)
    
    #Row minus1: Xlimits are the same as plus1
    Ym2=Y1-(3*Ydiff)
    
    #Row Plus2
    XpRP2=X2+(3*Xdiff);XmRP2=X1-(3*Xdiff)
    Yp3=Y1+(5*Ydiff)
    
    #Row minus2: Xlimits are the same as plus2
    Ym3=Y1-(5*Ydiff)
    
    #Row Plus3
    XpRP3=X2+(2*Xdiff);XmRP3=X1-(2*Xdiff)
    Yp4=Y1+(7*Ydiff)

    #Row minus3: Xlimits are the same as plus3
    Ym4=Y1-(7*Ydiff)
    #####END OF TOP LAYER BONDS
    
    #BOTTOM LAYER BOUNDS
    #[X1BOT,Y1BOT,Z1BOT]=removed_2_layers[metal][five-1]
    
    #bounds for middle Row BOTTOM layer
    XpRP0BOT=XpRP1; XmRP0BOT=XmRP1
    Yp1BOT=Y1BOT+Ydiff; Ym1BOT=Y1BOT-Ydiff
    
    #Row Plus1
    XpRP1BOT=XpRP2; XmRP1BOT=XmRP2
    Yp2BOT=Y1BOT+(3*Ydiff)
    
    #Row minus1: Xlimits are tha same as plus1
    Ym2BOT=Y1BOT-(3*Ydiff)
    
    #Row Plus2 
    XpRP2BOT=XpRP3; XmRP2BOT=XmRP3
    Yp3BOT=Y1BOT+(5*Ydiff)
    
    #Row minus2: Xlimits are tha same as plus2
    Ym3BOT=Y1BOT-(5*Ydiff)
    
    #Row Plus3 
    XpRP3BOT=X1BOT+(2*Xdiff); XmRP3BOT=X1BOT-(2*Xdiff)
    Yp4BOT=Y1BOT+(7*Ydiff)
    
    #Row minus3: Xlimits are tha same as plus3
    Ym4BOT=Y1BOT-(7*Ydiff)
    #END OF BOTTOM LAYER BOUBDS
    
    #print('Yp1=',Yp1,'Yp2=',Yp2,'Yp3=',Yp3,'YP4=',Yp4)
    #print('Ym1=',Ym1,'Ym2=',Ym2,'Ym3=',Ym3,'YP4=',Ym4)
    
    #print(XpRP0,XpRP1,XpRP2,XpRP3)
    #print(XmRP0,XmRP1,XmRP2,XmRP3)
    
    #Extract Cluster Atoms
    TopZero=[];TopPlus1=[];TopPlus2=[];TopPlus3=[];TopMinus1=[];TopMinus2=[];TopMinus3=[]
    BotZero=[];BotPlus1=[];BotPlus2=[];BotMinus1=[];BotMinus2=[];BotMinus3=[];BotPlus3=[]
    for atoms, coord in removed_2_layers.items():
        if atoms == metal:
            for i in range(0,len(coord)):
                Xcoord=removed_2_layers[atoms][i][0]
                Ycoord=removed_2_layers[atoms][i][1]
                Zcoord=removed_2_layers[atoms][i][2]
                #TOP LAYER
                if (Zcoord > limits[2]):
                    #Row Zero
                    if   ( Ym1 < Ycoord < Yp1 ) and ( XmRP0 < Xcoord < XpRP0 ):
                        TopZero.append(removed_2_layers[atoms][i])
                    #Row Plus1
                    elif ( Yp1 < Ycoord < Yp2 ) and ( XmRP1 < Xcoord < XpRP1 ):
                        TopPlus1.append(removed_2_layers[atoms][i])
                    #Row Plus2
                    elif ( Yp2 < Ycoord < Yp3 ) and ( XmRP2 < Xcoord < XpRP2 ):
                        TopPlus2.append(removed_2_layers[atoms][i])
                    #Row Plus3
                    elif ( Yp3 < Ycoord < Yp4 ) and ( XmRP3 < Xcoord < XpRP3 ):
                        TopPlus3.append(removed_2_layers[atoms][i])
                    #Row minus1
                    elif ( Ym2 < Ycoord < Ym1 ) and ( XmRP1 < Xcoord < XpRP1 ):
                        TopMinus1.append(removed_2_layers[atoms][i])
                    #Row minus2
                    elif ( Ym3 < Ycoord < Ym2 ) and ( XmRP2 < Xcoord < XpRP2 ):
                        TopMinus2.append(removed_2_layers[atoms][i])
                    #Row minus3
                    elif ( Ym4 < Ycoord < Ym3 ) and ( XmRP3 < Xcoord < XpRP3 ):
                        TopMinus3.append(removed_2_layers[atoms][i])
                #BOTTOM LAYER
                elif (Zcoord < limits[2]):
                    #Row Zero Bot
                    if ( Ym1BOT < Ycoord < Yp1BOT ) and ( XmRP0BOT < Xcoord < XpRP0BOT ):
                        BotZero.append(removed_2_layers[atoms][i])
                    #Row Plus1 Bot
                    if ( Yp1BOT < Ycoord < Yp2BOT) and ( XmRP1BOT < Xcoord < XpRP1BOT ):
                        BotPlus1.append(removed_2_layers[atoms][i])
                    #Row Plus2 Bot
                    elif ( Yp2BOT < Ycoord < Yp3BOT) and ( XmRP2BOT < Xcoord < XpRP2BOT ):
                        BotPlus2.append(removed_2_layers[atoms][i])
                    #Row minus1 Bot
                    elif ( Ym2BOT < Ycoord < Ym1BOT) and ( XmRP1BOT < Xcoord < XpRP1BOT ):
                        BotMinus1.append(removed_2_layers[atoms][i])
                    #Row minus2 Bot
                    elif ( Ym3BOT < Ycoord < Ym2BOT) and ( XmRP2BOT < Xcoord < XpRP2BOT ):
                        BotMinus2.append(removed_2_layers[atoms][i])
                    #Row minus3BOT, we need this row only if Y1BOT > Y1TOP
                    elif (Y1 < Y1BOT) and ( Ym4BOT < Ycoord < Ym3BOT) and ( XmRP3BOT < Xcoord < XpRP3BOT ):
                        BotMinus3.append(removed_2_layers[atoms][i])
                    #Row plus3BOT, we need this row only if Y1BOT < Y1TOP
                    elif (Y1 > Y1BOT) and ( Yp3BOT < Ycoord < Yp4BOT) and ( XmRP3BOT < Xcoord < XpRP3BOT ):
                        BotPlus3.append(removed_2_layers[atoms][i])
            

    #sort atoms in each row based on Xcoordinate and add them to cluster list
    cluster=[]       
    #BOTTOM LAYER
    if (Y1 < Y1BOT):
        BotMinus3.sort(key= lambda x: x[0]); 
        for i in BotMinus3: cluster.append(i)
        BotMinus2.sort(key= lambda x: x[0]); 
        for i in BotMinus2: cluster.append(i)
        BotMinus1.sort(key= lambda x: x[0]);
        for i in BotMinus1: cluster.append(i);
        BotZero.sort(key= lambda x: x[0]);
        for i in BotZero: cluster.append(i)
        BotPlus1.sort(key= lambda x: x[0]);
        for i in BotPlus1: cluster.append(i)
        BotPlus2.sort(key= lambda x: x[0]);
        for i in BotPlus2: cluster.append(i)
    elif (Y1 > Y1BOT):
        BotMinus2.sort(key= lambda x: x[0]);
        for i in BotMinus2: cluster.append(i)
        BotMinus1.sort(key= lambda x: x[0]);
        for i in BotMinus1: cluster.append(i);
        BotZero.sort(key= lambda x: x[0]);
        for i in BotZero: cluster.append(i)
        BotPlus1.sort(key= lambda x: x[0]);
        for i in BotPlus1: cluster.append(i)
        BotPlus2.sort(key= lambda x: x[0]);
        for i in BotPlus2: cluster.append(i)
        BotPlus3.sort(key= lambda x: x[0]);
        for i in BotPlus3:  cluster.append(i)
        
    #TOP LAYER
    TopMinus3.sort(key= lambda x: x[0]);
    for i in TopMinus3: cluster.append(i);
    TopMinus2.sort(key= lambda x: x[0]);
    for i in TopMinus2: cluster.append(i)
    TopMinus1.sort(key= lambda x: x[0]);
    for i in TopMinus1:  cluster.append(i);
    TopZero.sort(key= lambda x: x[0]);
    for i in TopZero: cluster.append(i)
    TopPlus1.sort(key= lambda x: x[0]);
    for i in TopPlus1: cluster.append(i);
    TopPlus2.sort(key= lambda x: x[0]);
    for i in TopPlus2: cluster.append(i)
    TopPlus3.sort(key= lambda x: x[0]);
    for i in TopPlus3: cluster.append(i)
        
    cluster_dict=removed_2_layers
    cluster_dict[metal]=cluster      
    
    tot_atoms=0 # number of atoms in cluster for header of xyz file
    for atoms, coord in cluster_dict.items():
        #print(np.shape(coord))
        tot_atoms += np.shape(coord)[0]
        
    
    #write xyz file
    #write header
    fwrite=open(OUT,'w')
    fwrite.write(str(tot_atoms) + '\n')
    fwrite.write('Cluster' + '\n')
    #write coordinates
    for key, values in cluster_dict.items():
        for i in range(0,len(values)):
            a=['{0:.15f}'.format(x) for x in cluster_dict[key][i]]
            fwrite.write(key+'        '+'      '.join(a)+'\n')
    fwrite.close()
