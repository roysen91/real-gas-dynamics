'''
Created on 16.02.2017

@author: roysonntag
'''
import matplotlib as mpl
import matplotlib.pyplot as plt
#mpl.use('TkAgg')
import cantera as ct
import numpy as np
import sys
sys.path.insert(0, 'lib')
import classes as cl
import combust as cmb
import station as stat
from constants import FARStoi, AirToFuelFlow,PressureConversion, FLowConversion,StandartPressure,StandartTemperature

def plot_combustion_temperature(t_in,p_in,w_in,comp,models=['RRD','CEA','CANTERA'],baseline='CANTERA',eqr_range=np.linspace(0,2,50)):
    fig = plt.figure()
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    legend1 = []
    legend2 = []

    if baseline == 'CANTERA':
        spec = ['O2','N2','CO2','AR','H2O']
        structure_cantera=''
        for s in spec:
            structure_cantera= structure_cantera+s+':'+str(comp.get_fraction(s))+','
        structure_cantera=structure_cantera[:-1]    
        # Get all of the Species objects
        species = {S.name: S for S in ct.Species.listFromFile('nasa.cti')}
        species2  = {S.name: S for S in ct.Species.listFromFile('gri30.xml')}
        # Create an IdealGas object including incomplete combustion species
        complete_species = [species2[S] for S in species2]
        complete_species.append(species['Jet-A(g)'])
        gas2 = ct.Solution(thermo='IdealGas', species=complete_species)#species.values())
        t_base = np.zeros(eqr_range.shape)
        #comb_air = 'O2:'+std_air.get_fraction('O2')
        for i,eqr in enumerate(eqr_range):
            gas2.TP = t_in, p_in
            gas2.set_equivalence_ratio(eqr, 'Jet-A(g)',structure_cantera)
            gas2.equilibrate('HP')
            
            t_base[i] = gas2.T
        # plot baseline values
        ax1.plot(eqr_range,t_base)
        legend1.append('Cantera')
    elif baseline == 'Online CEA':
        eqr_range_cea = []
        t_onl_cea=[]
        file1 = 'data/online_cea_combustion/largeChemSet_EQR_below1.dat'
        file2 = 'data/online_cea_combustion/largeChemSet_EQR_above1.dat'
        file1_handle = open(file1,'r')
        for line in file1_handle:
            if 'EQ.RATIO' in line:
                line=line.strip().split(',')
                eqr_range_cea.append(float(line[-1].split('=')[-1]))
            elif 'T, K' in line:
                line=line.strip().split(' ')
                t_onl_cea.append(float(line[-1]))   
        file1_handle.close()         
        file2_handle = open(file2,'r')
        for line in file2_handle:
            if 'EQ.RATIO' in line:
                line=line.strip().split(',')
                eqr_range_cea.append(float(line[-1].split('=')[-1]))
            elif 'T, K' in line:
                line=line.strip().split(' ')
                t_onl_cea.append(float(line[-1]))   
        file2_handle.close() 
        eqr_range = np.array(eqr_range_cea)
        t_base    = np.array(t_onl_cea)
        # plot baseline values
        ax1.plot(eqr_range,t_base)
        legend1.append('Online CEA') 
                                     
                                     

    # set up combustor
    comb = cmb.Combustor()
    comb.in_station = stat.Station(t_in,p_in,w_in,comp)
    comb.out_station = stat.Station(t_in,p_in,w_in,comp)
    comb.current_fuel.tf = 382.5948
    # get input fuel flow
    far_range = eqr_range*FARStoi
    wf_range = far_range*comb.in_station.w*AirToFuelFlow
    # get RRD values
    if 'RRD' in models:
        # set up combustor for RRD calculation
        comb.current_fluid = cl.RRDFluid()
        comb.current_fuel.burn_model = 'RRD'
        t_rrd = []
        # get t_out for every fuel flow
        for fuel_flow in np.nditer(wf_range):
            comb.wf = fuel_flow
            comb.heat_balance_calcs()
            t_rrd.append(comb.out_station.t)
        t_rrd = np.array(t_rrd)
        diff_t_rrd = (t_rrd-t_base)/t_base
        ax1.plot(eqr_range,t_rrd)
        ax2.plot(eqr_range,diff_t_rrd)
        legend1.append('RRD')
        legend2.append('RRD')
    
    # get CEA values
    if 'CEA' in models:
        # set up combustor for RRD calculation
        comb.current_fluid = cl.CeaFluid()
        comb.current_fuel.burn_model = 'CEA'
        t_cea = []
        # get t_out for every fuel flow
        for fuel_flow in np.nditer(wf_range):
            comb.wf = fuel_flow
            comb.heat_balance_calcs()
            t_cea.append(comb.out_station.t)
        t_cea = np.array(t_cea)
        diff_t_cea = (t_cea-t_base)/t_base
        ax1.plot(eqr_range,t_cea)
        ax2.plot(eqr_range,diff_t_cea)
        legend1.append('CEA')
        legend2.append('CEA')
        
    if 'CANTERA' in models:
        spec = ['O2','N2','CO2','AR','H2O']
        structure_cantera=''
        for s in spec:
            structure_cantera= structure_cantera+s+':'+str(comp.get_fraction(s))+','
        structure_cantera=structure_cantera[:-1]    
        # Get all of the Species objects
        species = {S.name: S for S in ct.Species.listFromFile('nasa.cti')}
        species2  = {S.name: S for S in ct.Species.listFromFile('gri30.xml')}
        # Create an IdealGas object including incomplete combustion species
        complete_species = [species2[S] for S in species2]
        complete_species.append(species['Jet-A(g)'])
        gas2 = ct.Solution(thermo='IdealGas', species=complete_species)#species.values())
        t_cantera = np.zeros(eqr_range.shape)
        #comb_air = 'O2:'+std_air.get_fraction('O2')
        for i,eqr in enumerate(eqr_range):
            gas2.TP = t_in, p_in
            gas2.set_equivalence_ratio(eqr, 'Jet-A(g)',structure_cantera)
            gas2.equilibrate('HP')
            
            t_cantera[i] = gas2.T
        diff_t_cantera = (t_cantera-t_base)/t_base
        # plot baseline values
        ax1.plot(eqr_range,t_cantera)
        legend1.append('Cantera')
        ax2.plot(eqr_range,diff_t_cantera)
        legend2.append('Cantera')
        
    ax1.set_xlabel('Equivalence ratio, $\phi$')
    ax1.set_ylabel('Temperature [K]')
    ax1.grid(True)
    ax1.legend(legend1)
    
    ax2.set_xlabel('Equivalence ratio, $\phi$')
    ax2.set_ylabel('delta T [%]')
    ax2.grid(True)
    ax2.legend(legend2)
    
    plt.show()
        
dummy_fluid = cl.CeaFluid()
std_air=cl.Composition('Air')
dummy_fluid.set_humidity(StandartTemperature+20, StandartPressure, std_air, RH=1)
plot_combustion_temperature(846,512*PressureConversion,57*FLowConversion,std_air,baseline='Online CEA')    