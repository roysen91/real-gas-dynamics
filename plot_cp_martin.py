'''
Created on 22 Feb 2017

@author: roysonntag
'''
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import classes as cl
import combust as cmb
import station as stat
from constants import StandartPressure,StandartTemperature,FLowConversion,PressureConversion,AirToFuelFlow

dummy_fluid = cl.CeaFluid()
comp=cl.Composition('Air')



humidity = [x/100 for x in range(0,100)]

# inlet conditions
t_in=846
p_in=512*PressureConversion
w_in=57*FLowConversion
# set up combustor
comb = cmb.Combustor()
comb.current_fuel.tf = 382.5948
comb.current_fluid = cl.CeaFluid()
comb.current_fuel.burn_model = 'CEA'


fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)
h2o=[]
o2=[]
co2=[]
n2=[]
ar=[]
h_in=[]
h_out=[]
for rh in humidity:
    dummy_fluid.set_humidity(StandartTemperature+20, StandartPressure, comp, RH=rh)
    # set stations with new comp 
    comb.in_station = stat.Station(t_in,p_in,w_in,comp)
    comb.out_station = stat.Station(t_in,p_in,w_in,comp)
    comb.wf = 0.033*comb.in_station.w*AirToFuelFlow
    comb.heat_balance_calcs()
    h2o.append(comp.get_fraction('H2O'))
    o2.append(comp.get_fraction('O2'))
    co2.append(comp.get_fraction('CO2'))
    n2.append(comp.get_fraction('N2'))
    ar.append(comp.get_fraction('AR'))
    h_in.append(comb.in_station.h)
    h_out.append(comb.out_station.h)
    # reset comp 
    comp=cl.Composition('Air')
ax1.plot(humidity,h2o,humidity,co2,humidity,o2,humidity,n2,humidity,ar)
ax1.legend(['H2O','CO2','O2','N2','AR'])
ax1.set_xlabel('relative humidity, [%]')
ax1.set_ylabel('mole fraction [%]')
ax1.grid(True)

ax2.plot(humidity,h_out,humidity,h_in)
ax2.legend(['Outlet Enthalpy','Inlet Enthalpy'])
ax1.set_xlabel('relative humidity, [%]')
ax2.set_ylabel('enthalpy [J/kg]')
ax2.grid(True)
    
plt.show()
    