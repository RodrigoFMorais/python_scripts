#########################################################################################
#       Code to calculate the generalized cordination number when subtracting atom     #
#       by Rodrigo Ferreira de Morais                                                  #
#       Date: 09/10/2014                                                               #
########################################################################################

#!/usr/bin/python

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
        if str(s)==t[i][0]:
            result=int(t[i][1])
            break
            
    return result
#----------------------------------------------------------------------------

########################### Main code #######################################

#---Definition of variables---

site=''
AtomRemoved=[]
AtomListSite=[]
tempAR=[]
AtomCoord=[]
AtomCoordC=[]
AtomCN=[]
vect=[]
AtomN=0.0
AtomT=''
CN=[]
ListaNB=[]
ListaNBCrude=[]
coordRemoved=0

param=3.3 #Pt-Pt maximum distance in A.

CNmaxtable=[
['s',37],
['h',39],
]
#-----------------------------

try:
   type=sys.argv[1] #Name of the Vasp CONTCAR or POSCAR used in the calculations
   type=sys.argv[2] #list of neighbors separeted by comma 

except:
    print 'The correct use of this script is: '
    print 'GCNR  file x,y,z_of_the removed_atom'
    sys.exit(2)

TempAR=sys.argv[2].split(",")


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

AtomT=linhas[5].split()[0]
AtomN=int(linhas[6].split()[0])

for i in range(9,9+AtomN):
    temp=[]
    temp.append(float(linhas[i].split()[0]))
    temp.append(float(linhas[i].split()[1]))
    temp.append(float(linhas[i].split()[2]))
    AtomCoord.append(temp)


# Direct to cartesian coordinates
for i in range(0,AtomN):
    temp=[]
    temp.append(AtomCoord[i][0]*vect[0][0]+AtomCoord[i][1]*vect[1][0]+AtomCoord[i][2]*vect[2][0])
    temp.append(AtomCoord[i][0]*vect[0][1]+AtomCoord[i][1]*vect[1][1]+AtomCoord[i][2]*vect[2][1])
    temp.append(AtomCoord[i][0]*vect[0][2]+AtomCoord[i][1]*vect[1][2]+AtomCoord[i][2]*vect[2][2])  
    AtomCoordC.append(temp)

temp=[]
# Direct to cartesian coordinates for the removed atom
temp.append(float(TempAR[0])*vect[0][0]+float(TempAR[1])*vect[1][0]+float(TempAR[2])*vect[2][0])
temp.append(float(TempAR[0])*vect[0][1]+float(TempAR[1])*vect[1][1]+float(TempAR[2])*vect[2][1])
temp.append(float(TempAR[0])*vect[0][2]+float(TempAR[1])*vect[1][2]+float(TempAR[2])*vect[2][2])
AtomRemoved.append(temp)

# Calculate the normal coordination number
for i in range(0,len(AtomCoordC)):
    temp=0
    for j in range(0,len(AtomCoordC)):
        if i<>j :
            dist=CalcDist(AtomCoordC[i],AtomCoordC[j])
            if dist<=param:
                temp=temp+1
    
    dist=CalcDist(AtomCoordC[i],AtomRemoved[0])
    if dist<=param:
        coordRemoved=coordRemoved+1
        AtomListSite.append(i)

    CN.append(temp)

if coordRemoved==8:
    site='s'
else:
    site='h'

for i in range(0,len(AtomListSite)):
    for j in range(0,len(AtomCoordC)):
        if int(AtomListSite[i])<>j :
            dist=CalcDist(AtomCoordC[int(AtomListSite[i])],AtomCoordC[j])
            if dist<=param:
                ListaNBCrude.append(j)

ListaNB=CleanList(ListaNBCrude,AtomListSite)

CNmax=FindMax(site,CNmaxtable)

temp=0.0
temp1=[]
for i in range(0,len(ListaNB)):
    temp=temp+CN[ListaNB[i]]
    temp1.append(CN[ListaNB[i]])
d={x:temp1.count(x) for x in temp1}
print d 
#print "The generalized coordination number is: ", temp/CNmax, CNmax, d, AtomListSite
#----------------------------------------------------------------------------
