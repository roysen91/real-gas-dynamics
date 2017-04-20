'''
Created on 13.02.2017

@author: roysonntag
'''
from classes import *
import numpy as np
from constants import *
import matplotlib as mpl
import matplotlib.pyplot as plt

my_fuel = Composition(comp_def='Jet-A1')
comp = Composition()
my_fluid = CeaFluid(eos_model='SOA')
my_fluid.eos.set_coeffs(comp)

t_range = np.linspace(StandartTemperature,2000,200)
p = StandartPressure


kappa_ideal=[]
kappa_soa=[]
kappa_real=[]
for t in t_range:
    kappa_ideal.append(my_fluid.tp2kappa(t, p, comp))
    kappa_soa.append(my_fluid.tp2kappa_eos(t, p, comp,mode='real'))
    kappa_real.append(my_fluid.tp2kappa_eos(t, p, comp))
    
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
#ax1.plot(pressure_ratio,eta_const,pressure_ratio,eta_mean_ideal,pressure_ratio,eta_mean_real)
ax1.plot(t_range,kappa_ideal,t_range,kappa_soa,t_range,kappa_real)
ax1.set_xlabel('Temperature [K]')
ax1.set_ylabel('kappa [1]')
ax1.grid(True)
ax1.legend(['ideal $\kappa$','SOA $\kappa$','real $\kappa$'])
plt.show()