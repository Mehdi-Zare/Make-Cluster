# I have two versions of making cluster code (Cluster.jpg) in this repository (Python and Shell script)

# Python:

## The calling function is "Make_Cluster_51"
### This gets the name of your input file, which is vasp poscar formatted, and give you a cluster in .xyz format. metal keyword is your metal name, mode can be eSMS or iSMS, and finally dX,dY,dZ are the distances between your atoms in X,Y and Z direction that you can get them easyly from your input file
### warning: this function can only create a specific shape cluster that consists of 51 metal atom for 111 and 0001 cleaved surfaces
### Example in eSMS:
### Make_Cluster_51(IN='CONTCAR', OUT='cluster.xyz',metal='Ni',X=4,Y=5,mode='esms',dX=2.5,dY=2.16,dZ=1.996), It is worth mentioning that in this mode you need to input the 5 middle atoms of cluster that you want. To get them, run this program with some arbitrary numbers and it will write for you file removed_2_layers.xyz. From that file you can get those numbers. (Look at image Cluster.jpg)

### Example in iSMS Mode:
### Make_Cluster_51(IN='CONTCAR', OUT='cluster.xyz',metal='Ni',X=7,Y=7,mode='isms',dX=2.5,dY=2.16,dZ=1.996)


# Shell:

### It is only written like eSMS mode of the python code


