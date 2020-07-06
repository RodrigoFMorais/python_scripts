###############################################################################
#       Code to average the Pt-O bond distance for Pt201/H2O698 interface     #
#       and  calculate the surface and volume of a Nanoparticle for each step #
#       by Rodrigo Ferreira de Morais                                         #
#       Last modification date: 19/02/2015                                    #
###############################################################################

#!/usr/bin/python

import sys
import math

###################### Functions #############################################
def CalcArea(a,b,c):
    #------------------------------------------------------------------------#
    #          Area =   1/2 |u x v|,                                         #
    #          where:                                                        #
    #                 u = (xb-xa,yb-ya,zb-za)                                #
    #                 v = (xc-xa,yc-ya,zc-za)                                #
    #------------------------------------------------------------------------#

    x=(c[2]-a[2])*(b[1]-a[1])-(b[2]-a[2])*(c[1]-a[1])
    y=(b[2]-a[2])*(c[0]-a[0])-(b[0]-a[0])*(c[2]-a[2])
    z=(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    area = 0.5 * math.sqrt((x)**2 + (y)**2 + (z)**2)
    return area

def CalcVol(a,b,c,CM):

    #------------------------------------------------------------------------#
    #          Volume =   1/6 (h.|u x v|),                                   #
    #          where:                                                        #
    #                 obs.: tran = origin translated to the center of mass   #
    #                 u = (xbtran-xatran,ybtran-yatran,zbtran-zatran)        #
    #                 v = (xctran-xatran,yctran-yatran,zctran-zatran)        #
    #                 h = (xatran-xCMtran,yatran-yCMtran,zatran-zCMtran)     #
    #------------------------------------------------------------------------#
    #print "------"
    #print a
    #print b
    #print c
    #print "------"

    atran=[]
    btran=[]
    ctran=[]

    p=0.0
    t=0.0
    vol=0.0

    u=[]
    v=[]

    atran.append(a[0]-CM[0])
    atran.append(a[1]-CM[1])
    atran.append(a[2]-CM[2])

    # btran=u
    btran.append(b[0]-CM[0])
    btran.append(b[1]-CM[1])
    btran.append(b[2]-CM[2])

    # ctran=v
    ctran.append(c[0]-CM[0])
    ctran.append(c[1]-CM[1])
    ctran.append(c[2]-CM[2])

    u.append(btran[0]-atran[0])
    u.append(btran[1]-atran[1])
    u.append(btran[2]-atran[2])

    v.append(ctran[0]-atran[0])
    v.append(ctran[1]-atran[1])
    v.append(ctran[2]-atran[2])

    p = atran[0]*u[1]*v[2] + atran[1]*u[2]*v[0] + atran[2]*u[0]*v[1]

    t = atran[2]*u[1]*v[0] + atran[0]*u[2]*v[1] + atran[1]*u[0]*v[2]

    vol = (1.0/6.0)*(p-t)

    if vol<0:
        vol=vol*-1

    return vol

#----------------------------------------------------------------------------

def truncated_octahedron(coordN,CSurf,coordC,CM,cor):
#truncated_octahedron(ListaNei[i],CSurf[i],AtomListC[i],CM[i],dPtPt)
    #---------------------------------------------------------------------------------#
    #   calculate the area and Volument for truncated_octahedron nanoparticle shape   #
    #---------------------------------------------------------------------------------#
    
    Result=[]
    #CSurf=[]
    InCnSurf=[]
    TSurfA=0.0
    TVol=0.0
    A111=0.0
    A100=0.0


    # Creat surface mesh of triangles according to the rules fo the truncated octahedron
    for i in range(0,len(coordN)):
            for j in range(0,len(coordN[i])-1):
                for k in range(j+1,len(coordN[i])):
                    dist=CalcDist(coordC[coordN[i][j]],coordC[coordN[i][k]])
		    if CSurf[i]==4 or CSurf[coordN[i][j]]==4 or CSurf[coordN[i][k]]==4 or (CSurf[i]==5 and CSurf[coordN[i][j]]==5 or CSurf[coordN[i][k]]==5):  
                        if dist<=4.5:
                                temp=[]
                                temp.append(CalcArea(coordC[i],coordC[coordN[i][k]],coordC[coordN[i][j]])/2)
                                temp.append(CalcVol(coordC[i],coordC[coordN[i][k]],coordC[coordN[i][j]],CM)/2)
                                Result.append(temp)
				A100=A100+temp[0]
                    else:
                        if dist<=3.3:
                            temp=[]
                            temp.append(CalcArea(coordC[i],coordC[coordN[i][k]],coordC[coordN[i][j]])/3)
                            temp.append(CalcVol(coordC[i],coordC[coordN[i][k]],coordC[coordN[i][j]],CM)/3)
                            Result.append(temp)
			    A111=A111+temp[0]

    # total surface area and volume:
    for i in range(0,len(Result)):
        TSurfA=TSurfA+Result[i][0]
        TVol=TVol+Result[i][1]

    temp=[]
    temp.append(TSurfA)
    temp.append(TVol)
    temp.append(A100)
    temp.append(A111)

    return temp
#----------------------------------------------------------------------------


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

site=0
vect=[]
AtomList=[]
AtomListC=[]
ListaNB=[]
ListaNBCrude=[]
CM=[]
CNSurf=[]
ListaNei=[]
RESULT=[]

dPtPt=3.4 #Pt-Pt maximum distance in A.
dPtO=2.6 #Pt-O maximum distance in A.

#-----------------------------

try:
   type=int(sys.argv[1]) #Number of XDATCAR

except:
    print 'The correct use of this script is: '
    print 'python Code.py Number_of_XDATCAR_file  '
    sys.exit(2)

# ----------- Get xyz
file=open("XDATCAR_run1","r")

linhas=file.readlines()

file.close()

    
for i in range(2,5):
    temp=[]
    temp.append(float(linhas[i].split()[0]))
    temp.append(float(linhas[i].split()[1]))
    temp.append(float(linhas[i].split()[2]))
    vect.append(temp)

for nf in range(1,type+1):
        AtomPt=0
        AtomO=0
        AtomH=0

	file=open("XDATCAR_run" + str(nf),"r")
	linhas=[]
	linhas=file.readlines()
	file.close()

	AtomPt=int(linhas[6].split()[0])
	AtomO=int(linhas[6].split()[1])
	AtomH=int(linhas[6].split()[2])+int(linhas[6].split()[3])
	TAtomp1=AtomPt+AtomO+AtomH+1
	inter=int((len(linhas)-7)/(TAtomp1))


	for li in range(0,inter):
		AtomPtCoord=[]
		AtomPtCoordC=[]
        	AtomOCoord=[]
        	AtomOCoordC=[]
        	AtomHCoord=[]
        	AtomHCoordC=[]
		DistPtO=[]
        	CNdPtO=[]
		CN=[]
		PtOmed=0.0
		Nsurf=0
		result=[]
		CNSurf=[]
        	ListaNei=[]
		CM=[]
		AtomListC=[]		

		tempx=0.0
		tempy=0.0
		tempz=0.0

		# Pt atoms and transform in XYZ
		for i in range(8+li*TAtomp1,8+AtomPt+li*TAtomp1):
    			temp=[]
    			temp.append(float(linhas[i].split()[0]))
    			temp.append(float(linhas[i].split()[1]))
    			temp.append(float(linhas[i].split()[2]))
    			AtomPtCoord.append(temp)

		for i in range(0,AtomPt):
    			temp=[]
    			temp.append(AtomPtCoord[i][0]*vect[0][0]+AtomPtCoord[i][1]*vect[1][0]+AtomPtCoord[i][2]*vect[2][0])
    			temp.append(AtomPtCoord[i][0]*vect[0][1]+AtomPtCoord[i][1]*vect[1][1]+AtomPtCoord[i][2]*vect[2][1])
    			temp.append(AtomPtCoord[i][0]*vect[0][2]+AtomPtCoord[i][1]*vect[1][2]+AtomPtCoord[i][2]*vect[2][2])  
    			AtomPtCoordC.append(temp)
			tempx=tempx+temp[0]
			tempy=tempy+temp[1]
			tempz=tempz+temp[2]
		CM.append(tempx/len(AtomPtCoordC))
		CM.append(tempy/len(AtomPtCoordC))
		CM.append(tempz/len(AtomPtCoordC))

		# O atoms and transform in XYZ
		for i in range(8+AtomPt+li*TAtomp1,8+AtomPt+AtomO+li*TAtomp1):
    			temp=[]
    			temp.append(float(linhas[i].split()[0]))
    			temp.append(float(linhas[i].split()[1]))
    			temp.append(float(linhas[i].split()[2]))
    			AtomOCoord.append(temp)

		for i in range(0,AtomO):
    			temp=[]
    			temp.append(AtomOCoord[i][0]*vect[0][0]+AtomOCoord[i][1]*vect[1][0]+AtomOCoord[i][2]*vect[2][0])
    			temp.append(AtomOCoord[i][0]*vect[0][1]+AtomOCoord[i][1]*vect[1][1]+AtomOCoord[i][2]*vect[2][1])
    			temp.append(AtomOCoord[i][0]*vect[0][2]+AtomOCoord[i][1]*vect[1][2]+AtomOCoord[i][2]*vect[2][2])
    			AtomOCoordC.append(temp)

		# H atoms and transform in XYZ
		for i in range(8+AtomPt+AtomO+li*TAtomp1,8+AtomPt+AtomO+AtomH+li*TAtomp1):
    			temp=[]
    			temp.append(float(linhas[i].split()[0]))
    			temp.append(float(linhas[i].split()[1]))
    			temp.append(float(linhas[i].split()[2]))
    			AtomHCoord.append(temp)

		for i in range(0,AtomH):
    			temp=[]
    			temp.append(AtomHCoord[i][0]*vect[0][0]+AtomHCoord[i][1]*vect[1][0]+AtomHCoord[i][2]*vect[2][0])
    			temp.append(AtomHCoord[i][0]*vect[0][1]+AtomHCoord[i][1]*vect[1][1]+AtomHCoord[i][2]*vect[2][1])
    			temp.append(AtomHCoord[i][0]*vect[0][2]+AtomHCoord[i][1]*vect[1][2]+AtomHCoord[i][2]*vect[2][2])
    			AtomHCoordC.append(temp)

		# calculate the CN for each Pt
		for i in range(0,len(AtomPtCoordC)):
    			temp=0
    			for j in range(0,len(AtomPtCoordC)):
        			if i<>j :
            				dist=CalcDist(AtomPtCoordC[i],AtomPtCoordC[j])
            				if dist<=dPtPt:
                				temp=temp+1
    			CN.append(temp)

		temp=[]
		for i in range(0,len(CN)):
			if CN[i]<10: 
				Nsurf=Nsurf+1
				temp.append(AtomPtCoord[i])
				AtomListC.append(AtomPtCoordC[i])
		AtomList.append(temp)

		temp4=[]
        	# calculate the CNsurf for each Pt
        	for i in range(0,len(AtomListC)):
                	temp=0
	        	temp3=[]
                	for j in range(0,len(AtomListC)):
                        	if i<>j :
                                	dist=CalcDist(AtomListC[i],AtomListC[j])
                                	if dist<=dPtPt:
                                        	temp=temp+1
						temp3.append(j)
			ListaNei.append(temp3)
			CNSurf.append(temp)


		# find the Pt-O bounds
		for i in range(0,len(AtomPtCoordC)):
    			for j in range(0,len(AtomOCoordC)):
        			dist=CalcDist(AtomPtCoordC[i],AtomOCoordC[j])
        			if dist<=dPtO:
    	   				DistPtO.append(dist)
	   				CNdPtO.append(CN[i])

		for i in range(0,len(DistPtO)):
    			PtOmed=PtOmed+DistPtO[i]
		PtOmed=PtOmed/len(DistPtO)

		result.append(truncated_octahedron(ListaNei,CNSurf,AtomListC,CM,dPtPt))

		temp=[]
		temp.append(len(DistPtO))
                temp.append(Nsurf)
                temp.append(PtOmed)
                temp.append(result[0][0])
                temp.append(result[0][1])
                temp.append(result[0][2])
                temp.append(result[0][3])
		RESULT.append(temp)

	
		#print li,len(DistPtO),Nsurf,PtOmed,result[li][0],result[li][1],result[li][2],result[li][3]

#Creat XDATCAR with the atoms of the surface.
Creatxdatcar(AtomList,vect,linhas[0])

for i in range(0,len(RESULT)):
	#     step,coverage,Nsurf,PtOmed,Area,volume,A111,A100
	print i,RESULT[i][0],RESULT[i][1],RESULT[i][2],RESULT[i][3],RESULT[i][4],RESULT[i][5],RESULT[i][6]


#----------------------------------------------------------------------------
