#!/bin/bash 

##################################################################################################################################
#				Author : Mehdi Zare on 3/4/2019                                                                  #
# 				Purpose: to make 51 atom cluster out of CONTCAR							 #
#
#
#################################################################################################################################

###  INSTRUCTION 
echo $'\n' This code need the CONTCAR file and create for you your cluster file containes 51 metal plus adsorbate $'\n'

# Check if CONTCAR exist
if [ ! -f "CONTCAR" ]; then echo $'\n' I have found CONTCAR file, please provide CONTCAR file  $'\n'; exit 1; fi

# Keep original CONTCAR
cp CONTCAR CONTCAR-original

sed  -i '/^[[:space:]]*$/,$d' CONTCAR

metalID=`head -6 CONTCAR | tail -1 | awk '{print $1}'`; metalnum=`head -7 CONTCAR | tail -1 | awk '{print $1}'`;
ads1ID=`head -6 CONTCAR | tail -1 | awk '{print $2}'` ; ads1num=`head -7 CONTCAR | tail -1 | awk '{print $2}'`;
ads2ID=`head -6 CONTCAR | tail -1 | awk '{print $3}'` ; ads2num=`head -7 CONTCAR | tail -1 | awk '{print $3}'`;
ads3ID=`head -6 CONTCAR | tail -1 | awk '{print $4}'` ; ads3num=`head -7 CONTCAR | tail -1 | awk '{print $4}'`;

echo $'\n' $metalID is $metalnum and $ads1ID is $ads1num and $ads2ID is $ads2num and and $ads3ID is $ads3num $'\n'

Adsorbnum=$(($ads1num+$ads2num+ads3num))

###First Part-Creating an expanded surface###
# Get the supercell size that you want
echo  $'\n' please enter multipliciy in X and Y direction that you are interested in, they must be integer numbers $'\n'
#####################################  FLAG ######################################
echo  $'\n' do you have these numbers, yes or no $'\n'
read answer
if [ "$answer" = "yes" ]
  then
      echo $'\n' alright, Lets continue $'\n' pelase enter X and Y in order $'\n'
elif [ "$answer" = "no" ]
  then
     echo $'\n' See you later, the program is going to exit! $'\n'
     exit 1
  else
     echo $'\n' please insert the correct word $'\n'
      exit 1
 fi
##################################  END OF FLAG ####################################
read xvalue
read yvalue


vasp-build-supercell-from-poscar -f CONTCAR -o XYZ -X $xvalue -Y $yvalue > expanded.xyz

# remove header
head -2 expanded.xyz > first-header
totatom=`head -1 expanded.xyz`   # total atoms in your supercell 
totmetals=$(($xvalue*$yvalue*$metalnum))
sed '1,2 d' expanded.xyz > No-header    # remove header of the xyz file



# Find the maximum Z coordinate
head -n$metalnum No-header > metals
Zmax=0
	line=1
	while [ "$line" -le "$metalnum"  ]; do
                Zcoord=`head -n$line metals | tail -1 | awk '{print $4}'`; echo Zcoord is $Zcoord
                if [ $(echo "$Zcoord > $Zmax" | bc) -eq 1 ]; then Zmax=$Zcoord; fi;
        line=$(($line+1))
        done
echo Zmax is $Zmax
top=$Zmax; bot=0

diff=`awk "BEGIN { print $top - $bot}"`;
interlayer=`awk "BEGIN { print $diff / 3 }"`; limit=`awk "BEGIN { print $interlayer / 2  }"`
 echo $'\n' the heigth difference of your bottom and top layer is $diff $'\n' and the avarage distance between layers is about $interlayer $'\n'
# set limits for differnet layers
limit1=`awk "BEGIN { print $bot + (0*$interlayer) + $limit }"`; limit2=`awk "BEGIN { print $bot + (1*$interlayer) + $limit  }"`;
limit3=`awk "BEGIN { print $bot + (2*$interlayer) + $limit  }"`
 echo $'\n' limits are $limit1 ,  $limit2 ,  $limit3 $'\n'

# tag atoms of the first two bottom layers with remove label
line=1
while [ "$line" -le "$totatom" ]; do
    Zcoord=`head -n$line No-header | tail -1 | awk '{print $4}'`; #### FALG #### echo Zcoord is $Zcoord
    if   [ $(echo "$Zcoord < $limit2" | bc) -eq 1  ]; then sed -i "$line s/$/ remove/" No-header ;fi
line=$(($line+1))
done




### remove all adsorbates and keeep adsorbate in the middile cell (13th cell)
varC=$ads1num       # number of carbon atoms of the adsorbate
varO=$ads2num       # number of oxygen atoms of the adsorbate
varH=$ads3num       # number of Hydrogen atoms of the adsorbate

cellnum=0
j=-1
trashcells=$((($xvalue*$yvalue)-1))
R=$((($xvalue*$yvalue)%2))                                # Check if the number of unitcells is an odd or even number
if [ "$R" -eq  0  ]; then theone=$((($xvalue*$yvalue)/2)); elif [ "$R" -ne  0  ]; then theone=$(((($xvalue*$yvalue)/2)+1)); fi    	# the middle cell number that we are interested in
while [ "$cellnum" -le "$trashcells"  ]                    # for cell #1 to cell #4 and cell number 6 to 9
   do cellnum=$(($cellnum+1))                             # the first cell is cell #1
      j=$(($j+1))                                         # the first index is zero
      if [ "$cellnum" -eq "$theone" ]; then continue; fi             # in order to skip cell #5, the cell that we are interested in
      Cidfirst=$(($totmetals+($j*$varC)+1)); Cidlast=$(($Cidfirst+$varC-1))                           # carbon atoms id starts from after all Pt atoms (64x9)=576
      Oidfirst=$(($totmetals+($xvalue*$yvalue*$varC)+($j*$varO)+1)); Oidlast=$(($Oidfirst+$varO-1))                # oxygen atoms id starts from after all Pt atoms (64x9)=576 and all carbon atoms (9x2)
      Hidfirst=$(($totmetals+($xvalue*$yvalue*$varC)+($xvalue*$yvalue*$varO)+($j*$varH)+1)); Hidlast=$(($Hidfirst+$varH-1)) # hydroge  atoms id starts from after all Pt atoms (64x9)=1280 and all carbon atoms (9x2) and all Oxygen(9x2)  
     echo cellnumber is $cellnum , Cidfirst is $Cidfirst and Cidlast is $Cidlast
     echo cellnumber is $cellnum , Oidfirst is $Oidfirst and Oidlast is $Oidlast
     echo cellnumber is $cellnum , Hidfirst is $Hidfirst and Hidlast is $Hidlast
   sed -i "$Cidfirst,$Cidlast s/$/ remove/" No-header                  # tag Carbon  corresponding atoms
   sed -i "$Oidfirst,$Oidlast s/$/ remove/" No-header                  # tag Oxygen corresponding atoms
   sed -i "$Hidfirst,$Hidlast s/$/ remove/" No-header                  # tag Hydrogen corresponding atoms
done

### end of remving adsorbates except for cell #10
sed  '/remove/d' No-header > Order.xyz            # remove all lines with remove tag, it gives us two top layer and adsorbate atoms of just one cell
Reallines=$((($totmetals/2)+$Adsorbnum)); Now=`wc -l Order.xyz | awk '{ print $1  }'`
#CHECK the number of lines in firle Order.xyz
if [ "$Reallines" -ne "$Now"  ]; then echo $'\n' The number of lines in file Order.xyz must be $Reallines but I found it $Now $'\n'; exit 1; fi


# Get atom id from user
# top layer
echo  $'\n' please enter two middle atoms of your desired cluster in x direction order from file Order.xyz $'\n'
#####################################  FLAG ######################################
echo  $'\n' do you have these numbers, yes or no $'\n'
read answer
if [ "$answer" = "yes" ]
  then
      echo $'\n' alright, Lets continue $'\n' pelase enter middle atom numbers $'\n'
elif [ "$answer" = "no" ]
  then
     echo $'\n' please open Order.xyz file in Vesta and choose the atom numbers $'\n'
     sed -i  "1 i ORDER"  Order.xyz;sed -i  "1 i $Now"  Order.xyz;
     exit 1
  else
     echo $'\n' please insert the correct word $'\n'
     sed -i  "1 i ORDER"  Order.xyz;sed -i  "1 i $Now"  Order.xyz;
      exit 1
 fi
##################################  END OF FLAG ####################################

read one
read two

# Get atom id from user
# top layer
echo  $'\n' please enter two middle atoms of your desired cluster in y direction order from file Order.xyz $'\n'
#####################################  FLAG ######################################
echo  $'\n' do you have these numbers, yes or no $'\n'
read answer
if [ "$answer" = "yes" ]
  then
      echo $'\n' alright, Lets continue $'\n' pelase enter middle atom numbers $'\n'
elif [ "$answer" = "no" ]
  then
     echo $'\n' please open Order.xyz file in Vesta and choose the atom numbers $'\n'
     sed -i  "1 i ORDER"  Order.xyz;sed -i  "1 i $Now"  Order.xyz;
     exit 1
  else
     echo $'\n' please insert the correct word $'\n'
     sed -i  "1 i ORDER"  Order.xyz;sed -i  "1 i $Now"  Order.xyz;
      exit 1
 fi
##################################  END OF FLAG ####################################

read three
read four

# FOR BOTTOM LAYER
# Get atom id from user
echo  $'\n' please enter the middle atom of bottom layer of your desired cluster between the two middle atoms of top layer from file Order.xyz $'\n'
#####################################  FLAG ######################################
echo  $'\n' do you have this number, yes or no $'\n'
read answer
if [ "$answer" = "yes" ]
  then
      echo $'\n' alright, Lets continue $'\n' pelase enter middle atom numbers $'\n'
elif [ "$answer" = "no" ]
  then
     echo $'\n' please open Order.xyz file in Vesta and choose the atom numbers $'\n'
     sed -i  "1 i ORDER"  Order.xyz;sed -i  "1 i $Now"  Order.xyz;
     exit 1
  else
     echo $'\n' please insert the correct word $'\n'
     sed -i  "1 i ORDER"  Order.xyz;sed -i  "1 i $Now"  Order.xyz;
      exit 1
 fi
##################################  END OF FLAG ####################################

read five

# TOP LAYER BOUNDS
X1=`head -n$one Order.xyz | tail -1 | awk '{print $2}'`; echo X1 is $X1; Y1=`head -n$one Order.xyz | tail -1 | awk '{print $3}'`; echo Y1 is $Y1
X2=`head -n$two Order.xyz | tail -1 | awk '{print $2}'`; echo X2 is $X2; Y2=`head -n$two Order.xyz | tail -1 | awk '{print $3}'`; echo Y2 is $Y2
X3=`head -n$three Order.xyz | tail -1 | awk '{print $2}'`; echo X3 is $X3;Y3=`head -n$three Order.xyz | tail -1 | awk '{print $3}'`; echo Y3 is $Y3
X4=`head -n$four Order.xyz | tail -1 | awk '{print $2}'`; Y4=`head -n$four Order.xyz | tail -1 | awk '{print $3}'`; echo Y4 is $Y4

Xdiff=`awk "BEGIN { print (($X2 - $X1) / 2) }"`; Ydiff=`awk "BEGIN { print (($Y4 - $Y1) / 2) }"`
Xcenter=`awk "BEGIN { print $Xdiff + $X1 }"`

# bounds for middle Row top layer
Yp1=`awk "BEGIN { print $Y1 + $Ydiff }"`; Ym1=`awk "BEGIN { print $Y1 - $Ydiff }"`;
XpRP0=`awk "BEGIN { print $X2 + (5 * $Xdiff) }"`; XmRP0=`awk "BEGIN { print $X1 - (5 * $Xdiff) }"`;

# Row plus1
Yp2=`awk "BEGIN { print $Y1 + (3 * $Ydiff) }"`;
XpRP1=`awk "BEGIN { print $X2 + (4 * $Xdiff) }"`; XmRP1=`awk "BEGIN { print $X1 - (4 * $Xdiff) }"`;

# Row minus1: Xlimits are the same as plus1
Ym2=`awk "BEGIN { print $Y1 - (3 * $Ydiff) }"`;

# Row plus2
Yp3=`awk "BEGIN { print $Y1 + (5 * $Ydiff) }"`;
XpRP2=`awk "BEGIN { print $X2 + (3 * $Xdiff) }"`; XmRP2=`awk "BEGIN { print $X1 - (3 * $Xdiff) }"`;

# Row minus2: Xlimits are the same as plus2
Ym3=`awk "BEGIN { print $Y1 - (5 * $Ydiff) }"`;

# Row plus3
Yp4=`awk "BEGIN { print $Y1 + (7 * $Ydiff) }"`;
XpRP3=`awk "BEGIN { print $X2 + (2 * $Xdiff) }"`; XmRP3=`awk "BEGIN { print $X1 - (2 * $Xdiff) }"`;

# Row minus3: Xlimits are the same as plus3
Ym4=`awk "BEGIN { print $Y1 - (7 * $Ydiff) }"`;
#### END OF TOP LAYER BOUNDS

# BOTTOM LAYER BOUNDS
X1BOT=`head -n$five Order.xyz | tail -1 | awk '{print $2}'`; echo X1BOT is $X1BOT; Y1BOT=`head -n$five Order.xyz | tail -1 | awk '{print $3}'`; echo Y1BOT is $Y1BOT

# bounds for middle Row BOTTOM layer
Yp1BOT=`awk "BEGIN { print $Y1BOT + $Ydiff }"`; Ym1BOT=`awk "BEGIN { print $Y1BOT - $Ydiff }"`;
XpRP0BOT=`awk "BEGIN { print $X1BOT + (5 * $Xdiff) }"`; XmRP0BOT=`awk "BEGIN { print $X1BOT - (5 * $Xdiff) }"`;

# Row plus1
Yp2BOT=`awk "BEGIN { print $Y1BOT + (3 * $Ydiff) }"`;
XpRP1BOT=`awk "BEGIN { print $X1BOT + (4 * $Xdiff) }"`; XmRP1BOT=`awk "BEGIN { print $X1BOT - (4 * $Xdiff) }"`;

# Row minus1: Xlimits are the same as plus1
Ym2BOT=`awk "BEGIN { print $Y1BOT - (3 * $Ydiff) }"`;

# Row plus2
Yp3BOT=`awk "BEGIN { print $Y1BOT + (5 * $Ydiff) }"`;
XpRP2BOT=`awk "BEGIN { print $X1BOT + (3 * $Xdiff) }"`; XmRP2BOT=`awk "BEGIN { print $X1BOT - (3 * $Xdiff) }"`;

# Row minus2: Xlimits are the same as plus2
Ym3BOT=`awk "BEGIN { print $Y1BOT - (5 * $Ydiff) }"`;

# Row plus3
Yp4BOT=`awk "BEGIN { print $Y1BOT + (7 * $Ydiff) }"`;
XpRP3BOT=`awk "BEGIN { print $X1BOT + (2 * $Xdiff) }"`; XmRP3BOT=`awk "BEGIN { print $X1BOT - (2 * $Xdiff) }"`;

# Row minus3: Xlimits are the same as plus3
Ym4BOT=`awk "BEGIN { print $Y1BOT - (7 * $Ydiff) }"`;
#### END OF BOTTOM LAYER BOUNDS

echo Yp1 is $Yp1, Yp2 is $Yp2, Yp3 is $Yp3 and Yp4 is $Yp4 
echo Ym1 is $Ym1, Ym2 is $Ym2, Ym3 is $Ym3 and Ym4 is $Ym4
echo XpRP0 is $XpRP0 and XpRP1 is $XpRP1 and XpRP2 is $XpRP2 and XpRP3 is $XpRP3
echo XpRP0 is $XmRP0 and XmRP1 is $XmRP1 and XmRP2 is $XmRP2 and XmRP3 is $XmRP3


if [ -f "cluster.xyz"  ]; then echo $'\n'  I found cluster.xyz file, this will ber overwritten by the new one $'\n'; rm cluster.xyz; fi
# loop for finding cluster atoms of top laye
line=1
while [ "$line" -le "$Reallines"  ]; do
   id=`head -n$line Order.xyz | tail -1 | awk '{print $1}'`;
   X=`head -n$line Order.xyz | tail -1 | awk '{print $2}'`;
   Y=`head -n$line Order.xyz | tail -1 | awk '{print $3}'`;
   Z=`head -n$line Order.xyz | tail -1 | awk '{print $4}'`;
   # check netal ID
   if [ "$id" == "$metalID"  ]; then
	# check to be in top layer
	if   [ $(echo "$Z > $limit3" | bc) -eq 1  ]; then
		# Row zero
		if   [ $(echo "$Y > $Ym1 && $Y < $Yp1 && $X > $XmRP0 && $X < $XpRP0" | bc) -eq 1  ]; then head -n$line Order.xyz | tail -1 >> TopZero.xyz;
		# Row plus1
		elif [ $(echo "$Y > $Yp1 && $Y < $Yp2 && $X > $XmRP1 && $X < $XpRP1" | bc) -eq 1  ]; then head -n$line Order.xyz | tail -1 >> TopPlus1.xyz;
		# Row plus2
		elif [ $(echo "$Y > $Yp2 && $Y < $Yp3 && $X > $XmRP2 && $X < $XpRP2" | bc) -eq 1  ]; then head -n$line Order.xyz | tail -1 >> TopPlus2.xyz;
		# Row plus3
		elif [ $(echo "$Y > $Yp3 && $Y < $Yp4 && $X > $XmRP3 && $X < $XpRP3" | bc) -eq 1  ]; then head -n$line Order.xyz | tail -1 >> TopPlus3.xyz;
		# Row minus1
		elif [ $(echo "$Y > $Ym2 && $Y < $Ym1 && $X > $XmRP1 && $X < $XpRP1" | bc) -eq 1  ]; then head -n$line Order.xyz | tail -1 >> TopMinus1.xyz;
		# Row minus2
		elif [ $(echo "$Y > $Ym3 && $Y < $Ym2 && $X > $XmRP2 && $X < $XpRP2" | bc) -eq 1  ]; then head -n$line Order.xyz | tail -1 >> TopMinus2.xyz;
		# Row minus3
		elif [ $(echo "$Y > $Ym4 && $Y < $Ym3 && $X > $XmRP3 && $X < $XpRP3" | bc) -eq 1  ]; then head -n$line Order.xyz | tail -1 >> TopMinus3.xyz;
		fi
	elif  [ $(echo "$Z < $limit3" | bc) -eq 1  ]; then
                # Row zeroBOT
                if   [ $(echo "$Y > $Ym1BOT && $Y < $Yp1BOT && $X > $XmRP0BOT && $X < $XpRP0BOT" | bc) -eq 1  ]; then head -n$line Order.xyz | tail -1 >> BotZero.xyz;
                # Row plus1BOT
                elif [ $(echo "$Y > $Yp1BOT && $Y < $Yp2BOT && $X > $XmRP1BOT && $X < $XpRP1BOT" | bc) -eq 1  ]; then head -n$line Order.xyz | tail -1 >> BotPlus1.xyz;
                # Row plus2BOT
                elif [ $(echo "$Y > $Yp2BOT && $Y < $Yp3BOT && $X > $XmRP2BOT && $X < $XpRP2BOT" | bc) -eq 1  ]; then head -n$line Order.xyz | tail -1 >> BotPlus2.xyz;
                # Row minus1BOT
                elif [ $(echo "$Y > $Ym2BOT && $Y < $Ym1BOT && $X > $XmRP1BOT && $X < $XpRP1BOT" | bc) -eq 1  ]; then head -n$line Order.xyz | tail -1 >> BotMinus1.xyz;
                # Row minus2BOT
                elif [ $(echo "$Y > $Ym3BOT && $Y < $Ym2BOT && $X > $XmRP2BOT && $X < $XpRP2BOT" | bc) -eq 1  ]; then head -n$line Order.xyz | tail -1 >> BotMinus2.xyz;
                # Row minus3BOT, we need this row only if Y1BOT > Y1TOP
                elif [ $(echo "$Y1BOT > $Y1 && $Y > $Ym4BOT && $Y < $Ym3BOT && $X > $XmRP3BOT && $X < $XpRP3BOT" | bc) -eq 1  ]; then head -n$line Order.xyz | tail -1 >> BotMinus3.xyz;
                # Row plus3BOT, we need this row only if Y1BOT < Y1TOP
                elif [ $(echo "$Y1BOT < $Y1 && $Y > $Yp3BOT && $Y < $Yp4BOT && $X > $XmRP3BOT && $X < $XpRP3BOT" | bc) -eq 1  ]; then head -n$line Order.xyz | tail -1 >> BotPlus3.xyz;	
		fi		

	fi
   fi
line=$(($line+1))
done

# put atoms in order
wc -l Bot* Top* 
for row in Bot* Top* ; do total=`wc -l $row | awk '{ print $1  }'`
  # based on the second field whic is X coordinate
  sort -k2n -o $row $row 
done
# if we have BotMinus3.xyz it should come first, if not BotPlus3.xyz  will go to the end of the bottom layer
if [ -f "BotMinus3.xyz" ]; then 
  cat BotMinus3.xyz BotMinus2.xyz BotMinus1.xyz BotZero.xyz BotPlus1.xyz BotPlus2.xyz TopMinus3.xyz TopMinus2.xyz TopMinus1.xyz TopZero.xyz TopPlus1.xyz TopPlus2.xyz TopPlus3.xyz > cluster.xyz
  rm BotMinus3.xyz BotMinus2.xyz BotMinus1.xyz BotZero.xyz BotPlus1.xyz BotPlus2.xyz TopMinus3.xyz TopMinus2.xyz TopMinus1.xyz TopZero.xyz TopPlus1.xyz TopPlus2.xyz TopPlus3.xyz
else
  cat BotMinus2.xyz BotMinus1.xyz BotZero.xyz BotPlus1.xyz BotPlus2.xyz BotPlus3.xyz TopMinus3.xyz TopMinus2.xyz TopMinus1.xyz TopZero.xyz TopPlus1.xyz TopPlus2.xyz TopPlus3.xyz > cluster.xyz
  rm BotMinus2.xyz BotMinus1.xyz BotZero.xyz BotPlus1.xyz BotPlus2.xyz BotPlus3.xyz TopMinus3.xyz TopMinus2.xyz TopMinus1.xyz TopZero.xyz TopPlus1.xyz TopPlus2.xyz TopPlus3.xyz
fi 

 # Add adsorbstes
tail -n$Adsorbnum Order.xyz  >> cluster.xyz

Reallines=$((51+$Adsorbnum)); Now=`wc -l cluster.xyz | awk '{ print $1  }'`
#CHECK the number of lines in firle cluster.xyz

if [ "$Reallines" -ne "$Now"  ]; then echo $'\n' The number of lines in file cluster.xyz must be $Reallines but I found it $Now $'\n'; exit 1; fi

sed -i  "1 i cluster"  cluster.xyz;sed -i  "1 i $Now"  cluster.xyz;

rm first-header metals No-header 
