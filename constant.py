####File with constant
import xml.etree.ElementTree as ET
root_node = ET.parse('vasprun.xml').getroot()

constantdict = {}
for tag in root_node.findall('parameters/separator'):
        key = tag.get('name')
        if key == 'electronic':
            for value in tag.findall('i'):
                key = value.get('name')
                constantdict[key] = value.text
