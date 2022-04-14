# import xml.dom.minidom

# xmldoc = xml.dom.minidom.parse('vasprun.xml')
# eigenvalues = xmldoc.getElementsByTagName('eigenvalues')

# print("%d eigenvalues" % eigenvalues.length)

# for item in eigenvalues:
#     band = item.getElementsByTagName('r')

# bandstructure = []      
# for i in band:
#     # print(i.firstChild.nodeValue)
#     bandstructure.append(i.firstChild.nodeValue)

# print(bandstructure[1])