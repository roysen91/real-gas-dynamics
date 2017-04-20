'''
Created on 30 Jan 2017

@author: u566064
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
import cantera as ct

#read in data from files
FileNames=['data/online_cea_combustion/reducedChemSet_EQR_below1.dat',
           'data/online_cea_combustion/reducedChemSet_EQR_above1.dat',
           'data/online_cea_combustion/largeChemSet_EQR_below1.dat',
           'data/online_cea_combustion/largeChemSet_EQR_above1.dat']

T = []
CP = []
EQR = []
    
for i,inFile in enumerate(FileNames):
    if i==2:
        Tred=T
        CPred=CP
        EQRred=EQR
        
        T = []
        CP = []
        EQR = []
    inputFile = open(inFile,'r')

    for line in inputFile:
        if 'T, K' in line:
            line = line.split()
            T.append(float(line[2]))
        if 'Cp, KJ/(KG)(K)' in line:
            line = line.split()
            CP.append(float(line[2]))
        if 'EQ.RATIO' in line:
            line = line.split(',')
            EQR.append(float(line[-1].split('=')[1]))

    inputFile.close()

################# plot data

# reduced set of species using my cea 
dummy_fluid = cl.CeaFluid()
std_air=cl.Composition('Air')
dummy_fluid.set_humidity(StandartTemperature+20, StandartPressure, std_air, RH=1)

comb = cmb.Combustor()
comb.current_fluid = cl.CeaFluid()
comb.in_station = stat.Station(846,512*PressureConversion,57*FLowConversion,std_air)
comb.out_station = stat.Station(846,512*PressureConversion,57*FLowConversion,std_air)
comb.current_fuel.tf = 382.5948

eqr_vec = np.linspace(0.25, 2, 20)
far_vec = eqr_vec*FARStoi
wf_vec = far_vec*comb.in_station.w*AirToFuelFlow

t_cea = []

for fuel_flow in np.nditer(wf_vec):
    comb.wf = fuel_flow
    comb.heat_balance_calcs()
    t_cea.append(comb.out_station.t)
    
comb.current_fluid = cl.RRDFluid()
comb.current_fuel.burn_model = 'RRD'

t_rrd = []

for fuel_flow in np.nditer(wf_vec):
    comb.wf = fuel_flow
    comb.heat_balance_calcs()
    t_rrd.append(comb.out_station.t)
    
t_cea = np.array(t_cea)
t_rrd = np.array(t_rrd)




# reduced set of species using online cea
Tred=np.array(Tred)
CPred=np.array(CPred)
EQRred=np.array(EQRred)

# full set of species
T=np.array(T)
CP=np.array(CP)
EQR=np.array(EQR)


diff_CP=(CPred-CP)/CP
diff_T=(Tred-T)/T
diff_T_cea = (np.interp(EQR,eqr_vec,t_cea)-T)/T
diff_T_rrd = (np.interp(EQR,eqr_vec,t_rrd)-T)/T

######## CANTERA calculation ############
# Get all of the Species objects defined in the GRI 3.0 mechanism
species = {S.name: S for S in ct.Species.listFromFile('nasa.cti')}
# Create an IdealGas object including incomplete combustion species
complete_species = [species[S] for S in ('Jet-A(g)','O2','N2','CO2','H2O','Ar')]
gas2 = ct.Solution(thermo='IdealGas', species=complete_species)#species.values())
T_incomplete = np.zeros(EQR.shape)
for i in range(len(EQR)):
    gas2.TP = 846, 512*PressureConversion
    gas2.set_equivalence_ratio(EQR[i], 'Jet-A(g)', 'O2:0.19775868763, N2:0.7371626995, AR:0.008841156551, CO2:0.0003011563203 , H2O: 0.05593629993')
    gas2.equilibrate('HP')
    T_incomplete[i] = gas2.T

fig = plt.figure()

ax = fig.add_subplot(1,2,1)
ax.plot(EQRred,Tred)
ax.plot(EQR,T)  
ax.plot(eqr_vec,t_cea)  
ax.plot(eqr_vec,t_rrd)  
ax.plot(EQR,T_incomplete)  

plt.xlim(0.25, 1.5) 
plt.xlabel('Equivalence ratio, $\phi$')
plt.ylabel('Temperature [K]')
plt.grid(True)
plt.legend(['Onl. CEA RCS','Onl. CEA full CS', 'My CEA','RRD','Cantera'])

ax =  fig.add_subplot(1,2,2)
ax.plot(EQRred,diff_T)
ax.plot(EQR,diff_T_cea)
ax.plot(EQR,diff_T_rrd)
plt.grid(True)

plt.xlabel('Equivalence ratio, $\phi$')
plt.ylabel('delta T')

plt.xlim(0.25, 1.5) 
plt.legend(['Onl. CEA RCS', 'My CEA','RRD'])

plt.show()

