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
RRFluid = Fluid()
print(RRFluid.tp2h(300,101300,air))


##air.compress(p1,p2,T1,T2)