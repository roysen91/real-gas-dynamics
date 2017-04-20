'''
Created on 10.02.2017

@author: roysonntag
'''
import matplotlib as mpl
from constants import FARStoi, AirToFuelFlow, PressureConversion, FLowConversion,SpecHeatConversion,\
    StandartTemperature, StandartPressure
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import classes as cl
import combust as cmb
import station as stat


# to be implemented in Model
def get_ref_cp(t,p,comp):
    '''
    Reads JANAF tables for cp  for given temperature range. Interpolates
    if needed. Cp of mixture is calculated using volume weighted average
    INPUT:  list     t     Temperature
            float    p     Pressure
            obj      comp  Composition
    OUTPUT: list     cp    Spec. Heat Constant from JANAF tables
    '''
    # list holding dict for each species
    species_list = []
    
    for spec,frac in comp.structure.items():
        file = open('data/janaf_tables/janaf_'+spec.name.lower()+'_g_1bar.txt','r')
        data_array = np.loadtxt(file,skiprows=3)
        # only store temperature and cp values, name and fraction
        interp_arry = np.zeros([len(t),2])
        interp_arry[:,0] = t
        # interpolate read data for given temp range
        interp_arry[:,1] = np.interp(t,data_array[:,0],data_array[:,1])
        species_dict = {'data':interp_arry, 'name':spec.name, 'fraction':frac, 'mol_wgt': spec.mol_wgt}
        # append species dict to list
        species_list.append(species_dict)
        file.close()
    
    cp_array = mix_ideal(species_list)
    # convert from [J/K*mol] to [J/K*kg]
    cp_array[:,1] /= comp.mol_wgt#*comp.mol_wgt
    
    # return interpolated for given temperature range
    return cp_array[:,1]

def mix_ideal(species_list):
    '''
    calculates ideal mixture properties
    INPUT:  list         data        list of dict's for every species (keys: data,name,fraction)
    OUTPUT: np.array     mix_data    janaf data for ideal mixture
    '''    
    
    # init zero array of max size
    mix_data = np.zeros([len(species_list[0]['data'][:,0]),2])
    # iterate through data dict and fractions
    mix_data[:,0] = species_list[0]['data'][:,0]
    for sp in species_list:
        #mix_data[:,1] += sp['fraction']*sp['data'][:,1]*sp['mol_wgt']
        mix_data[:,1] += sp['fraction']*sp['data'][:,1]
    return mix_data 
        

################################## CP Plot
t_range = np.linspace(300,2000,100)

p=1e5
fl_cea = cl.CeaFluid()
fl_rrd = cl.RRDFluid()

std_air= cl.Composition('Air')

fl_cea.set_humidity(StandartTemperature+20, StandartPressure, std_air, RH=1)
#std_air.set_WGR(0.0368527)


cea_cp = []
rrd_cp = []

for t in t_range:
    cea_cp.append(fl_cea.tp2cp(t, p, std_air))
    rrd_cp.append(fl_rrd.tp2cp(t, p, std_air))
janaf_cp = get_ref_cp(t_range, p, std_air)
cea_cp = np.array(cea_cp)
rrd_cp = np.array(rrd_cp)

diff_cp_cea = (cea_cp-janaf_cp)/janaf_cp
diff_cp_rrd = (rrd_cp-janaf_cp)/janaf_cp

fig_cp = plt.figure()

ax = fig_cp.add_subplot(1,2,1)
ax.plot(t_range,janaf_cp)
ax.plot(t_range,cea_cp)
ax.plot(t_range,rrd_cp)

plt.legend(['JANAF','CEA', 'RRD'])

ax =  fig_cp.add_subplot(1,2,2)
ax.plot(t_range,diff_cp_cea)
ax.plot(t_range,diff_cp_rrd)
plt.legend(['CEA', 'RRD'])

plt.show()

