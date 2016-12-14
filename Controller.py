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


print('T_out:',RRFluid.isentropic_compression(288.15,101325,500000,air))
print('T_out:',BFluid.isentropic_compression(288.15,101325,500000,air))
print('T_out:',RFluid.isentropic_compression(288.15,101325,500000,air))

print('T_out:',RRFluid.polytropic_compression(0.9,288.15,101325,500000,air))
print('T_out:',BFluid.polytropic_compression(0.9,288.15,101325,500000,air))
print('T_out:',RFluid.polytropic_compression(0.9,288.15,101325,500000,air))
#print('RRD_h',RFluid.tp2h(300,101325,air))
#print('RRD_s',RFluid.tp2s(300,101325,air))
#print('RRD_cp',RFluid.tp2cp(300,101325,air))

gas = Composition()


##air.compress(p1,p2,T1,T2)