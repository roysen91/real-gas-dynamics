'''
This is the main plot file. It imports all necessary libs and calls various plotting functions
'''
import sys
sys.path.insert(0, 'lib')
import plot_functions as plotfcn
import classes as cl
import numpy as np
from constants import PressureConversion, FLowConversion,StandartPressure,StandartTemperature


# init dummy fluid class and standart air as composition 
dummy_fluid = cl.CeaFluid()
std_air 	=cl.Composition('Air')
# set relative humidity of air to 1
dummy_fluid.set_humidity(StandartTemperature+20, StandartPressure, std_air, RH=1)

#plotfcn.plot_combustion_temperature(846,512*PressureConversion,57*FLowConversion,std_air,baseline='Online CEA')

#plotfcn.plot_cp()

pressure_ratio = np.linspace(1,100,400)
#plotfcn.plot_work_press_ratio(pressure_ratio,StandartTemperature)
#plotfcn.plot_eff_press_ratio(pressure_ratio,StandartTemperature)
plotfcn.plot_burn_temp_comparison(512,846,56,t_fuel=382.5948)