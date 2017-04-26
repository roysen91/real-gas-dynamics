'''
Created on 13.02.2017

@author: roysonntag
'''
import sys
sys.path.insert(0, 'lib')
from classes import *
import numpy as np
from constants import *
import matplotlib as mpl
import matplotlib.pyplot as plt

my_fuel = Composition(comp_def='Jet-A1')
comp = Composition()
#my_fluid = CeaFluid(eos_model='SOA')
#my_fluid.eos.set_coeffs(comp)

eqr_range = np.linspace(0.1,1,10)
p = StandartPressure

#get CEA tool generated data
t_onl_cea=[]
file1 = 'data/CEA_tool/cea_full_rh=1.dat'
file1_handle = open(file1,'r')
for line in file1_handle:
    line=line.strip().split(' ')
    t_onl_cea.append(float(line[-1].split('=')[-1]))
file1_handle.close()   

# calc t_out using own routine      

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
#ax1.plot(pressure_ratio,eta_const,pressure_ratio,eta_mean_ideal,pressure_ratio,eta_mean_real)
ax1.plot(t_range,kappa_ideal,t_range,kappa_soa,t_range,kappa_real)
ax1.set_xlabel('Temperature [K]')
ax1.set_ylabel('kappa [1]')
ax1.grid(True)
ax1.legend(['ideal $\kappa$','SOA $\kappa$','real $\kappa$'])
plt.show()