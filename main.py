import xml.etree.ElementTree as ET

root_node = ET.parse('vasprun.xml').getroot()
# for child in root_node:
#     print(child.tag, child.attrib)

for tag in root_node.findall('kpoints/varray'):
    value = tag.get('name')
    if value is not None: print(value)

kpointslist = []

for tag in root_node.findall('kpoints/varray'):
    value = tag.get('name')
    if value == 'kpointlist':
        for i in tag.findall('v'):
            kpointslist.append(i.text)
            print(i.text)

bandlist = []

for num, tag in enumerate(root_node.findall('calculation/eigenvalues/array/set/set/set'), start=1):
    value = tag.get('comment')
    if value == 'kpoint {}'.format(num):
        for i in tag.findall('r'):
            bandlist.append(i.text)
            print(i.text)