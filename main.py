from functools import reduce
from pickle import TRUE
from re import S
import xml.etree.ElementTree as ET
import func
import constant
import numpy as np


HARTREE = 0.036749308136649
BOHR = 0.5291


root_node = ET.parse('vasprun303015.xml').getroot()
# for child in root_node:
#     print(child.tag, child.attrib)


#parse structure
basis = []
rec_basis = []

for tag in root_node.findall('structure/crystal/varray'):
    value = tag.get('name')
    if value == 'basis':
        for v, i in enumerate(tag.findall('v'), start=0):
            basis.append(list(map(float ,i.text.split())))
        continue
    elif value == 'rec_basis':
        for i in tag.findall('v'):
            rec_basis.append(list(map(float ,i.text.split())))
        break

# BASIS, LATTIECE CONSTANT 
basis = np.asanyarray(basis)
rec_basis = np.asanyarray(rec_basis)
#LAT CONST
a=np.sqrt(np.sum(basis[0,:]**2))
b=np.sqrt(np.sum(basis[1,:]**2))
c=np.sqrt(np.sum(basis[2,:]**2))
#REC LAT CONST
b1 = np.sqrt(np.sum(rec_basis[0,:]**2))
b2 = np.sqrt(np.sum(rec_basis[1,:]**2))
b3 = np.sqrt(np.sum(rec_basis[2,:]**2))
#Angl of lattiece
alpha = np.arccos(np.dot(basis[1,:],basis[2,:])/b/c)*360/2/np.pi
beta  = np.arccos(np.dot(basis[0,:],basis[2,:])/a/c)*360/2/np.pi
gamma = np.arccos(np.dot(basis[0,:],basis[1,:])/a/b)*360/2/np.pi

#BLOCK WHICH CONVERT TO RELEATIVE VALUES
reduced_a = a
reduced_b = b / a
reduced_c = c / a

reduced_lattiece_constant = np.array([reduced_a/BOHR, reduced_b, reduced_c])
reduced_basis = np.divide(basis, (reduced_a))
reduced_rec_basis = np.multiply(rec_basis, (2 * np.pi * reduced_a))

#K-POINTS parsing
kpointslist = []
temp = []
for tag in root_node.findall('kpoints/varray'):
    value = tag.get('name')
    if value == 'kpointlist':
        for i in tag.findall('v'):
            kpointslist.append(i.text[6:]) #append without space
            temp.append(list(map(float ,i.text.split())))


count_kx = 0
count_ky = 0
count_kz = 0
Nz = 0

for i in temp:
    for j, v in enumerate(i, start=1):
        if j == 1 and v != 0:
            count_kx  += 1
            v += v / b1 
        elif j == 2 and v != 0:
            count_ky  += 1
            v += v / b2 
        elif j == 3 and v != 0:
            count_kz  += 1

print(int(b1/0.03333333*3))
print(int(b2/0.03333333*3))
print(int((b3*(1/0.5))+0.5))

# BAND pasing
bandlist = []

for num, tag in enumerate(root_node.findall('calculation/eigenvalues/array/set/set/set'), start=1):
    value = tag.get('comment')
    if value == 'kpoint {}'.format(num):
        for i in tag.findall('r'):
            # bandlist.append(i.text[:-9]) #append without ~1 and space
            # print(float(i.text[0:-11])*HARTREE)
            bandlist.append(str(float(i.text[0:-11])*HARTREE)+'    ')
            

#spliting one list to len(list)/NBANDS shape
mergedlist = list(func.list_split(bandlist, int(constant.constantdict['NBANDS'])))

#write in file
with open('band.BND', 'w') as f:
    np.savetxt(f, (reduced_lattiece_constant), fmt="%.6f",newline=' ')
    np.savetxt(f, (reduced_basis), fmt="%.6f")
    np.savetxt(f, (reduced_rec_basis),fmt="%.6f")
    for i, (kpoint, band) in enumerate(zip(kpointslist, mergedlist),start=1):
         f.writelines(kpoint + '        /' + str(i) + '\n')
         f.writelines("".join(band)+'\n') # conver list to str

