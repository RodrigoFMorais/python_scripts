#######################################################################################
#   Code to average calculate the Velocity-Velocity autocorrelation function          #
#       by Rodrigo Ferreira de Morais                                                 #
#       Last modification date: 04/03/2015                                            #
#                                                                                     #
#       Zm_bar= (1/(M-m)) * Sum{n=0,n=M-m-1} < V[t(n+m)].V[t(n)]>/< V[t(0)].V[t(0)]>  #
#       where :                                                                       #
#              M = total number of steps that will be take in to account              #
#              Comp. Mat. Sci. 2 221-232 (1994)                                       #
#######################################################################################

#!/usr/bin/python

import sys
import math

###################### Functions #############################################

def Creatxdatcar(AtomList,vectors,Name):

    File=open("NSurfXDATCAR","w\n")

    File.write(Name)
    File.write("  1\n")
    File.write("  %.6f %.6f %.6f\n"% (vectors[0][0],vectors[0][1],vectors[0][2]))
    File.write("  %.6f %.6f %.6f\n"% (vectors[1][0],vectors[1][1],vectors[1][2]))
    File.write("  %.6f %.6f %.6f\n"% (vectors[2][0],vectors[2][1],vectors[2][2]))
    File.write("  Pt \n")
    File.write("  %i\n"%len(AtomList[0]))
    for i in range(0,len(AtomList)):
    	File.write("Direct configuration=     %i\n"%(i+1))
	for j in range(0,len(AtomList[i])):
        	tmp=[]
        	tmp.append(AtomList[i][j][0])
        	tmp.append(AtomList[i][j][1])
        	tmp.append(AtomList[i][j][2])
        	File.write("  %.8f %.8f %.8f \n"% (tmp[0],tmp[1],tmp[2]))

    return

def CalcDist(a,b):
    dist = math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)
    return dist

def CleanList(a,b):
    # ListaNBCrude,AtomListSite
  
    list2 = []    
    [list2.append(i) for i in a if not i in list2]

    for i in range(0,len(b)):
        list2.remove(int(b[i]))

    return list2

def FindFunc(a,b,n):
    # value_to_look,list_to_look, position_already_in
    result=1
    for i in range(0,len(b)):
        if i<>n:
            if b[i]==a:
                result=2
                break
    return result


def FindMax(s,t):
    result=0

    for i in range(0,len(t)):
        if s==t[i][0]:
            result=int(t[i][1])
            break
            
    return result
#----------------------------------------------------------------------------

########################### Main code #######################################

#---Definition of variables---

TAtom=0
Stpoint=1000
TSteps=0
vect=[]
Zm=[]
#-----------------------------

try:
   type=sys.argv[1] #Number of VDAT

except:
    print 'The correct use of this script is: '
    print 'python CalculateZn.py VDAT  '
    sys.exit(2)

# ----------- Get xyz
file=open(type,"r")

linhas=file.readlines()

file.close()

    
for i in range(2,5):
    temp=[]
    temp.append(float(linhas[i].split()[0]))
    temp.append(float(linhas[i].split()[1]))
    temp.append(float(linhas[i].split()[2]))
    vect.append(temp)

TAtom=int(linhas[6].split()[0])+int(linhas[6].split()[1])+int(linhas[6].split()[2])+int(linhas[6].split()[3])

TSteps=(len(linhas)-8)/(TAtom+5)

for i in range(Stpoint,TSteps-1):
	temp1=0.0
	for j in range(0,TSteps-i-1):
		for k in range(0,TAtom):
			for ind in range(0,3):
				temp1=temp1+float(linhas[8+(TAtom+5)*(i+j)+k].split()[ind])+float(linhas[8+(TAtom+5)*j+k].split()[ind]) 
	if (i==Stpoint) :
		Zm.append(temp1/(TSteps-i))
	else :
		
		Zm.append(temp1/(TSteps-i/Zm[0]))	


for i in range(0,len(Zm)):
	print i,Zm[0] 

#----------------------------------------------------------------------------
