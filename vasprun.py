from cmath import sqrt
from functools import reduce
from pickle import TRUE
from re import S
from tkinter import N
import xml.etree.ElementTree as ET
import func

import numpy as np



HARTREE = 0.036749308136649
BOHR = 0.5291

class Vasprun:
    """parser for vasprun.xml"""

    def __init__(self, vaspfile):
        self.root_node = ET.parse(vaspfile).getroot()
        
        self.__basis = []
        self.__recbasis = []
        self.__eigenvalues = []
        self.__kpoints = []

    @property
    def a(self):
        return self.__a

    @property
    def b(self):
        return self.__b

    @property
    def c(self):
        return self.__c

    @property
    def b1(self):
        return self.__b1

    @property
    def b2(self):
        return self.__b2

    @property
    def b3(self):
        return self.__b3

    @property
    def alpha(self):
        return self.__alpha

    @property
    def beta(self):
        return self.__beta

    @property
    def gamma(self):
        return self.__gamma


    def get_structure(self):
        """
        Parsing of structural constants 
        """
        
        for tag in self.root_node.findall('structure'):
            if tag.get('name') == 'finalpos':
                value = tag.findall('crystal/varray')
                for i in value:
                    if i.get('name') == 'basis':
                        for j in i:
                            self.__basis.append(list(map(float ,j.text.split())))
                    elif i.get('name') == 'rec_basis':
                        for j in i:
                            self.__recbasis.append(list(map(float ,j.text.split())))  

        # BASIS, LATTIECE CONSTANT 
        self.__basis = np.asanyarray(self.__basis)
        self.__recbasis = np.asanyarray(self.__recbasis)    

        #LAT CONST
        self.__a=np.sqrt(np.sum(self.__basis[0,:]**2))
        self.__b=np.sqrt(np.sum(self.__basis[1,:]**2))
        self.__c=np.sqrt(np.sum(self.__basis[2,:]**2))
        #REC LAT CONST
        self.__b1 = np.sqrt(np.sum(self.recbasis[0,:]**2))
        self.__b2 = np.sqrt(np.sum(self.recbasis[1,:]**2))
        self.__b3 = np.sqrt(np.sum(self.recbasis[2,:]**2))
        #Angl of lattiece
        self.__alpha = np.arccos(np.dot(self.__basis[1,:],self.__basis[2,:])/self.__b/self.__c)*360/2/np.pi
        self.__beta  = np.arccos(np.dot(self.__basis[0,:],self.__basis[2,:])/self.__a/self.__c)*360/2/np.pi
        self.__gamma = np.arccos(np.dot(self.__basis[0,:],self.__basis[1,:])/self.__a/self.__b)*360/2/np.pi


    @property
    def eigenvalues(self):
        return self.__eigenvalues
    
    @property
    def basis(self):
        return self.__basis

    @property
    def recbasis(self):
        return self.__recbasis


    def get_eigenvalues(self):
        """
        Parsing eigenvalues for band structure in list format
        """

        for num, tag in enumerate(self.root_node.findall('calculation/eigenvalues/array/set/set/set'), start=1):
            value = tag.get('comment')
            if value == 'kpoint {}'.format(num):
                for i in tag.findall('r'):
                    self.__eigenvalues.append(float(i.text[0:-11]))

    @property
    def kpoints(self):
        return self.__kpoints

    def get_kpoints(self):
        """
        Parsing kpoints
        """

        self.__kpoints = []

        for tag in self.root_node.findall('kpoints/varray'):
            value = tag.get('name')
            if value == 'kpointlist':
                for item in tag.findall('v'):
                    self.__kpoints.append(item.text.split())
                    for i in self.__kpoints:
                        for j,item in enumerate(i):
                            i[j] = float(item)


    def get_electronicparam(self):
        """_parsing of electronic parameters for bandtransform function_

        Returns:
            _dict_: _string_
        """        

        electronicparam = {}

        for tag in self.root_node.findall('parameters/separator'):
            key = tag.get('name')
            if key == 'electronic':
                for value in tag.findall('i'):
                    key = value.get('name')
                    electronicparam[key] = value.text

        return electronicparam

    @property
    def efermi(self):
        return self.__efermi

    @efermi.setter
    def efermi(self, value):
        print("Be careful!\nYou're changing the efermi\n")
        self.__efermi = value

    def get_efermi(self):
        self.__efermi = self.root_node.findall('calculation/dos/i')
        self.__efermi = [float(i.text) for i in self.__efermi]
        self.__efermi = self.__efermi[0]
        
    def get_subdivision(self):
        self.subdivision = []
        for tag in self.root_node.findall('kpoints/generation/v'):
            key = tag.get('name')
            if key == 'genvec1':
                self.subdivision.append(list(map(float, tag.text.split())))
            elif key == 'genvec2':
                self.subdivision.append(list(map(float, tag.text.split())))
            elif key == 'genvec3':
                self.subdivision.append(list(map(float, tag.text.split())))
        print(self.subdivision)

        self.subdivision = np.asanyarray(self.subdivision)
        self.h1=np.sqrt(np.sum(self.subdivision[0,:]**2))
        self.h2=np.sqrt(np.sum(self.subdivision[1,:]**2))
        self.h3=np.sqrt(np.sum(self.subdivision[2,:]**2))

        self.N1 = int((self.b1*3)/self.h1)
        self.N2 = int((self.b2*3)/self.h2)
        self.N3 = int((self.b3*3)/self.h3)
                


        
class Vasp2Igor(Vasprun):
    """
    Conevert Vasp out to Igor input
    """
    def __init__(self, vaspfile):
        super().__init__(vaspfile)
        super().get_structure()
        super().get_efermi()
        super().get_eigenvalues()
        super().get_kpoints()
        super().get_electronicparam()
        super().get_subdivision()
        


    def reducebasis(self):
        self.__reduced_a = self.__a
        self.__reduced_b = self.__b / self.__a
        self.__reduced_c = self.__c / self.__a

        self.__reduced_lattiece_constant = np.array([self.__reduced_a/BOHR, self.__reduced_b, self.__reduced_c])
        self.__reduced_basis = np.divide(self.__basis, (self.__reduced_a))
        self.__reduced_rec_basis = np.multiply(self.__rec_basis, (2 * np.pi * self.__reduced_a))


    def bandstransform(self):
        """Reading and replace eigevalues to eig*=HARTEE"""

        transformedbands = []

        for i in super().eigenvalues:
            transformedbands.append(i*HARTREE)

        transformedbands = list(func.list_split(transformedbands, int(self.get_electronicparam()["NBANDS"])))

        return transformedbands


    def transformkpoints(self):

        transformedkpoints = []

        for kpoints in super().kpoints:
            for kpoint in kpoints:
                transformedkpoints.append(kpoint*(2/1.73250))
        transformedkpoints = list(func.list_split(transformedkpoints, 3))

        return transformedkpoints

    
    def writefile(self, name):
        with open(name, 'w') as f:
            f.write(str(self.efermi))
            f.write()
            for i, (kpoint, band) in enumerate(zip(self.transformkpoints(), self.bandstransform()),start=1):
                # f.writelines(str(kpoint)[1:-1]+ '        /' + str(i) + '\n')
                # f.writelines(str(band)[1:-1]+'\n')
                f.write('  '.join(map(str, kpoint))+ '        /' + str(i) + '\n')
                f.write('  '.join(map(str, band))+'\n')



if __name__ == "__main__":
    # v = Vasp2Igor('vasprunrandom.xml')
    # v.get_kpoints()
    # v.get_eigenvalues()
    # v.get_structure()
    # v.writefile('tets2')
    v = Vasp2Igor('ignore/vasprun14127.xml')
    
    print(v.N1)
    print(v.N2)
    print(v.N3)
    # print(v.alpha, '\n', v.beta, '\n', v.gamma, '\n', v.a, '\n', v.b, '\n', v.c, '\n', v.b1, v.b2, v.b3)
    
    # N1 =14
    # N2 =12
    # N3 =7

    # print(v.b1/N1*3)
    # print(v.b2/N2*3)
    # print(v.b3/N3*3)
    
    
    # print(int((v.b1*3)/0.0726))
    # print(int((v.b2*3)/0.08333333))
    
