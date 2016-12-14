'''
Created on 17.11.2016

@author: roysonntag
'''
from classes import *
#import matplotlib.pyplot as plt
#plt.plot([1,2,3,4])
#plt.ylabel('some numbers')
#plt.show()


# init FLUID object
air = Composition()
RRFluid = CeaFluid()
BFluid = BuckerFluid()
RFluid = RRDFluid()

print('CEA_h',RRFluid.tp2h(300,101325,air))
print('CEA_s',RRFluid.tp2s(300,101325,air))
print('CEA_cp',RRFluid.tp2cp(300,101325,air))

print('Bucker_h',BFluid.tp2h(300,101325,air))
print('Bucker_s',BFluid.tp2s(300,101325,air))
print('Bucker_cp',BFluid.tp2cp(300,101325,air))

print('RRD_h',RFluid.tp2h(300,101325,air))
print('RRD_s',RFluid.tp2s(300,101325,air))
print('RRD_cp',RFluid.tp2cp(300,101325,air))

gas = Composition()


##air.compress(p1,p2,T1,T2)