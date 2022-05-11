from re import S
import xml.etree.ElementTree as ET
import func
import constant
import numpy as np

HARTREE = 0.036749308136649
BOHR = 0.5291

root_node = ET.parse('vasprun_dos.xml').getroot()
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

basis = np.asanyarray(basis)
rec_basis = np.asanyarray(rec_basis)
a=np.sqrt(np.sum(basis[0,:]**2))
b=np.sqrt(np.sum(basis[1,:]**2))
c=np.sqrt(np.sum(basis[2,:]**2))
print(a,b,c)

alpha = np.arccos(np.dot(basis[1,:],basis[2,:])/b/c)*360/2/np.pi
beta  = np.arccos(np.dot(basis[0,:],basis[2,:])/a/c)*360/2/np.pi
gamma = np.arccos(np.dot(basis[0,:],basis[1,:])/a/b)*360/2/np.pi
print(alpha,beta,gamma) 

basis = np.divide(basis, BOHR)
rec_basis = np.divide(rec_basis, BOHR)



#K-POINTS parsing
kpointslist = []

for tag in root_node.findall('kpoints/varray'):
    value = tag.get('name')
    if value == 'kpointlist':
        for i in tag.findall('v'):
            kpointslist.append(i.text[6:]) #append without space
            # print(i.text)

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
# mergedlist = list(func.list_split(bandlist, int(constant.constantdict['NBANDS'])))

#write in file
# with open('test.BND', 'w') as f:
#     for i in basis:
#         f.writelines(i + '\n')
#     for i in rec_basis:
#         f.writelines(i + '\n')
#     for i, (kpoint, band) in enumerate(zip(kpointslist, mergedlist),start=1):
#         f.writelines(kpoint + '        /' + str(i) + '\n')
#         f.writelines("".join(band)+'\n') # conver list to str

