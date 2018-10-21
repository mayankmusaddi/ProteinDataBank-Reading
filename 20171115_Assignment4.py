import numpy as n
import math
import sys

fname = input()
f = open(fname,'r')
oname=fname[:-4]+"_output.txt"
sys.stdout = open(oname,'w')
line = f.readline() #to ignore HEADER
line = f.readline() #to ignore HEADER

#TITLE printing------------------------------------
while line[0:6] == "TITLE ":
	print(line[10:80].partition("  ")[0],end = "")
	line = f.readline()
print()
#--------------------------------------------------

#CHAIN storage-------------------------------------
chain=[]
while line[0:6] != "COMPND":
	line = f.readline()
while line[0:6] == "COMPND":
	if line[11:17] == "CHAIN:":
		chain.append(line[18:80].partition(";")[0])
	line = f.readline()
chain=sorted(chain)
i=0
for key in chain:
	if len(key) > 1:
		chain[i]="("+key+")"
	i+=1
#--------------------------------------------------

#SEQRES storage------------------------------------
amin = {}
length = 0
while line[0:6] != "SEQRES":
	line = f.readline()
while line[0:6] == "SEQRES":
	seq = (line[19:70].split(" "))
	for key in seq :
		if key != '':
			length+=1
			if key in amin.keys():
				amin[key]+=1
			else:
				amin[key]=1
	line = f.readline()
print("LENGTH\t",length)
print("CHAIN\t",len(chain),"\t",','.join(str(i) for i in chain))

for key in sorted(amin.keys()):
	print(key,"\t",amin[key]/length)

if "UNK" in amin.keys():
	print("UNKNOWN\t",amin["UNK"])
else:
	print("UNKNOWN\t0")
#--------------------------------------------------

#HETNAM printing-----------------------------------
ligands=[]
while line[0:6] != "HETNAM":
	if line[0:6] == "CRYST1":
		break
	line = f.readline()
while line[0:6] == "HETNAM":
	ligands.append(line[11:14])
	line = f.readline()
ligands=sorted(list(set(ligands)))
if "H2O" in ligands:
	ligands.remove("H2O")
if len(ligands) != 0:
	print("LIGANDS\t",','.join(str(i) for i in ligands))
#-------------------------------------------------
#Function for Calculating Dihedral Angle
def calang(l1,l2,l3,l4):
	a=n.array(l1)
	b=n.array(l2)
	c=n.array(l3)
	d=n.array(l4)
	v1=a-b
	v2=b-c
	v3=c-d
	l1=n.cross(v1,v2)
	l2=n.cross(v2,v3)
	n1=l1/n.sqrt(n.dot(l1,l1))
	n2=l2/n.sqrt(n.dot(l2,l2))
	u1=v2/n.sqrt(n.dot(v2,v2))
	u2=n.cross(u1,n2)
	rad=-math.atan2(n.dot(n1,u2),n.dot(n1,n2))
	angle=n.degrees(rad)
	return (round(angle,3))

#ATOM angle---------------------------------------
while line[0:6] != "ATOM  ":
	line = f.readline()
cid=''
k=0
flag=0
while line[0:6] != "END   ":
	if line[0:6] == "ATOM  ":
		if line[21] != cid:
			cid = line[21]
			l1=[]
			l2=[]
			l3=[]
			if flag == 1:
				print("NA\tNA")
			print("CHAIN-",cid)
			flag=1
			print("NA",end="\t")
			k=1
		atom=line[13:16].split(" ")[0]
		x=float(line[30:38])
		y=float(line[38:46])
		z=float(line[46:54])
		pt=[x,y,z]
		if atom == "C":
			if len(l1) != 0:
				print(calang(l1,l2,l3,pt),end="\t")
				k+=1
			l1 = pt
		if atom == "N":
			if len(l2) != 0:
				print(calang(l2,l3,l1,pt),end="\t")
				k+=1
			l2 = pt
		if atom == "CA":
			if len(l3) != 0:
				print(calang(l3,l1,l2,pt),end="\t")
				k+=1
			l3 = pt
		if k%3 == 0:
			print()
	line = f.readline()

print("NA\tNA")
# ------------------------------------------------

f.close()
sys.stdout.close()