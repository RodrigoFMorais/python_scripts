###############################################################################
#       Code to find the surface atoms bounded to the O2                      #
#       by Rodrigo Ferreira de Morais                                         #
#       Date: 29/01/2015                                                      #
###############################################################################

#!/usr/bin/python

import os
import sys
import math

###################### Functions #############################################
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

site=0
AtomListSite=[]
AtomCoord=[]
AtomCoordC=[]
AtomCN=[]
vect=[]
AtomN=0.0
AtomT=''
CN=[]
ListaNB=[]
ListaNBCrude=[]

param=2.8 #Pt-Pt maximum distance in A.

#-----------------------------

try:
   type=sys.argv[1] #Name of the Vasp CONTCAR or POSCAR used in the calculations

except:
    print 'The correct use of this script is: '
    print 'FindN  file'
    print 'Ex.:'
    print 'FindN CONTCAR_run3'
    sys.exit(2)

# ----------- Get xyz
file=open(sys.argv[1],"r")

linhas=file.readlines()

file.close()

for i in range(2,5):
    temp=[]
    temp.append(float(linhas[i].split()[0]))
    temp.append(float(linhas[i].split()[1]))
    temp.append(float(linhas[i].split()[2]))
    vect.append(temp)

#AtomT=linhas[5].split()[0]
AtomN=0
for i in range(0,len(linhas[5].split())):
    AtomN=AtomN+int(linhas[5].split()[i])

for i in range(8,8+AtomN):
    temp=[]
    temp.append(float(linhas[i].split()[0]))
    temp.append(float(linhas[i].split()[1]))
    temp.append(float(linhas[i].split()[2]))
    AtomCoord.append(temp)

for i in range(0,AtomN):
    temp=[]
    temp.append(AtomCoord[i][0]*vect[0][0]+AtomCoord[i][1]*vect[1][0]+AtomCoord[i][2]*vect[2][0])
    temp.append(AtomCoord[i][0]*vect[0][1]+AtomCoord[i][1]*vect[1][1]+AtomCoord[i][2]*vect[2][1])
    temp.append(AtomCoord[i][0]*vect[0][2]+AtomCoord[i][1]*vect[1][2]+AtomCoord[i][2]*vect[2][2])  
    AtomCoordC.append(temp)

# Find the number of sites.
for j in range(0,len(AtomCoordC)-2):
    	dist=CalcDist(AtomCoordC[len(AtomCoordC)-2],AtomCoordC[j])
    	if dist<=param:
        	AtomListSite.append(j)

for j in range(0,len(AtomCoordC)-2):
        dist=CalcDist(AtomCoordC[len(AtomCoordC)-1],AtomCoordC[j])
        if dist<=param:
                AtomListSite.append(j)

#for i in range(0,len(AtomListSite)):
#	print AtomListSite[i]+1, 
#print len(AtomCoordC)-1, len(AtomCoordC)

for i in range(0,len(AtomListSite)):
 	os.system('mkdir Pt' + str(AtomListSite[i]+1))
#	os.system('cp ./DOS' + str(AtomListSite[i]+1) + '  ./Pt' + str(AtomListSite[i]+1) + '/.')
        os.system('/Users/rferre01/bin/python_codes/script_plot_dos DOS' + str(AtomListSite[i]+1) + ' ' + str(AtomListSite[i]+1) + ' Pt')
        os.system('mv *.eps ./Pt' + str(AtomListSite[i]+1) + '/') 

os.system('mkdir O1')
os.system('/Users/rferre01/bin/python_codes/script_plot_dos DOS' + str(len(AtomCoordC)-1) + ' ' + str(len(AtomCoordC)-1) + ' O')
os.system('mv *.eps ./O1/')  

os.system('mkdir O2')
os.system('/Users/rferre01/bin/python_codes/script_plot_dos DOS' + str(len(AtomCoordC)) + ' ' + str(len(AtomCoordC)) + ' O')
os.system('mv *.eps ./O2/') 



#os.system('cp *.eps ./GRAPHIC/')


#----------------------------------------------------------------------------
