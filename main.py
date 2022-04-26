import xml.etree.ElementTree as ET
import func
import constant
import numpy as np

root_node = ET.parse('vasprun.xml').getroot()

for child in root_node:
    print(child.tag, child.attrib)

#K-POINTS parsing
kpointslist = []

for tag in root_node.findall('kpoints/varray'):
    value = tag.get('name')
    if value == 'kpointlist':
        for i in tag.findall('v'):
            kpointslist.append(i.text)
            # print(i.text)

# BAND pasing
bandlist = []

for num, tag in enumerate(root_node.findall('calculation/eigenvalues/array/set/set/set'), start=1):
    value = tag.get('comment')
    if value == 'kpoint {}'.format(num):
        for i in tag.findall('r'):
            bandlist.append(i.text[:-9]) #append without ~1 
            # print(i.text)

#spliting one list to len(list)/NBANDS shape
mergedlist = list(func.list_split(bandlist, int(constant.constantdict['NBANDS'])))
with open('test.BND', 'w') as f:
    for i, (kpoint, band) in enumerate(zip(kpointslist, mergedlist),start=1):
        f.writelines(kpoint + '        /' + str(i) + '\n')
        f.writelines("".join(band)+'\n') # conver list to str

