# This file is inteded to demonstrate the workflow of my real gas library.
# It will calculate a combustion process for standart air with relative humidity 
# of 100% at DTAMB = 20K
# Baseline will be the result from the CEA desktop app from file "data/CEA_tool/cea_full_rh=1.dat"
# There will be two subplots, the first holding absoulute values and the second the deviation
# from the CEA desktop app results


# Import all nessecary libraries
# sys let's you import folders into working directory
import sys
sys.path.insert(0, 'lib')
import classes as cl 	# import gas properties
import combust as cmb
import station as stat
import numpy as np 		# numerical python for vector calcs
from constants import *		# import all constants
# import matplotlib environment for plotting
import matplotlib as mpl
import matplotlib.pyplot as plt

####### INPUTS #######
t_in = 846 # temperature [K]
p_in = 512 # pressure [psi]
w_in = 56 # air flow [lb/s]
dtamb = 20
rh = 1

# init dummy fluid class and standart air as composition 
dummy_fluid = cl.CeaFluid()
comp		= cl.Composition('Air')
# set relative humidity of air to 1
dummy_fluid.set_humidity(StandartTemperature+dtamb, StandartPressure,comp, RH=rh)

####### CALCULATIONS #######
# baseline for the plot will be online CEA results for rh=1
# initiate empty numpy arrays
eqr_range 	= np.array([])
t_range		= np.array([])

# declare file handle
file1='data/CEA_tool/cea_full_rh=1.dat'
# open file handle
with open(file1,'r') as file1_handle:
    # skip first row
    next(file1_handle)
    # read line by line
    for line in file1_handle:
    	# seperate each value
        line=line.strip().split('  ')
        # convert eqr and t to float and append to array
        eqr_range = np.append(eqr_range,float(line[0]))
        t_range   = np.append(t_range,float(line[1]))
    # close file handle
    file1_handle.close()

# set up combustor
comb 				= cmb.Combustor()
# station at inlet of combustor
comb.in_station 	= stat.Station(t_in,p_in,w_in,comp)
# station at outlet of combustor
comb.out_station 	= stat.Station(t_in,p_in,w_in,comp)

# get input fuel flow
far_range 			= eqr_range*FARStoi
wf_range 			= far_range*comb.in_station.w*AirToFuelFlow

# get RRD values
# set up combustor for RRD calculation
comb.current_fluid 	= cl.RRDFluid()
comb.current_fuel.burn_model = 'RRD'
t_rrd = np.array([])
# get t_out for every fuel flow
for fuel_flow in np.nditer(wf_range):
    comb.wf 	= fuel_flow
    comb.heat_balance_calcs()
    t_rrd 		= np.append(t_rrd,comb.out_station.t)
# get deviation from baseline value
diff_t_rrd 		= (t_rrd-t_range)/t_range*100

# get CEA values
# set up combustor for RRD calculation
comb.current_fluid = cl.CeaFluid()
comb.current_fuel.burn_model = 'CEA'
t_cea = np.array([])
# get t_out for every fuel flow
for fuel_flow in np.nditer(wf_range):
    comb.wf 	= fuel_flow
    comb.heat_balance_calcs()
    t_cea 		= np.append(t_cea,comb.out_station.t)
# get deviation from baseline value
diff_t_cea = (t_cea-t_range)/t_range*100


####### PLOTTING #######
fig = plt.figure() 				# invoke new figure 
ax1 = fig.add_subplot(1,2,1)	# first subplot
ax2 = fig.add_subplot(1,2,2)	# second subplot

# declare legends
legend1 = ['Online CEA','RRD','my CEA']
legend2 = ['RRD','my CEA']

# plot absolute values in first axis
ax1.plot(eqr_range,t_range,color='black')
ax1.plot(eqr_range,t_rrd)
ax1.plot(eqr_range,t_cea)
ax1.set_xlabel('Equivalence ratio, $\phi$')
ax1.set_ylabel('Temperature [K]')
ax1.grid(True)
ax1.legend(legend1)

# plot deviation from baseline in second axis
ax2.plot(eqr_range,diff_t_rrd)
ax2.plot(eqr_range,diff_t_cea)
ax2.set_xlabel('Equivalence ratio, $\phi$')
ax2.set_ylabel('delta T [%]')
ax2.set_ylim([0,5])
ax2.grid(True)
ax2.legend(legend2)

# show resulting figure
plt.show()