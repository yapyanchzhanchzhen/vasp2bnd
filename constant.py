####File with constant
import xml.etree.ElementTree as ET

root_node = ET.parse('vasprun303015.xml').getroot()

#PARSE ELECTRINIC CONSTANT
constantdict = {}
for tag in root_node.findall('parameters/separator'):
        key = tag.get('name')
        if key == 'electronic':
            for value in tag.findall('i'):
                key = value.get('name')
                constantdict[key] = value.text


# class VaspVariebles():
#     def __init__(self, PREC = None, ENMAX = None, ENAUG = None, EDIFF = None, IALGO = None, IWAVPR = None, NBANDS = None, NELECT = None, TURBO = None, IRESTART = None, \
#          NREBOOT = None, NMIN = None, EREF = None):
#         self.prec = PREC
#         self.enmax = ENMAX
#         self.enaug = ENAUG
#         self.ediff = EDIFF
#         self.ialgo = IALGO
#         self.iwavpr = IWAVPR
#         self.nbands = NBANDS
#         self.nelect = NELECT
#         self.turbo = TURBO
#         self.irestart = IRESTART
#         self.nreboot = NREBOOT
#         self.nmin = NMIN
#         self.eref = EREF

#     def parser(root_node):

#         constantdict = {}

#         for tag in root_node.findall('parameters/separator'):
#             key = tag.get('name')
#             if key == 'electronic':
#                 for value in tag.findall('i'):
#                     key = value.get('name')
#                     constantdict[key] = value.text


