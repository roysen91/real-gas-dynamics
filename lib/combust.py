'''
Created on 06.02.2017

@author: roysonntag
'''

import station as st
import classes as cl
import fuel
from constants import *

class Combustor():
    
    
    def __init__(self):
        self.eta_b = 1 
        self.wf = 0                 # fuel flow
        self.fgrc = 0               # corrected fuel to gas ratio
        self.in_station = st.Station() 
        self.out_station = st.Station() 
        self.current_fuel = fuel.Fuel() 
        self.current_fluid = cl.CeaFluid()
        
    def compute(self):
        pass
    
    def heat_balance_calcs(self):
        self.fgrc = self.wf/(self.in_station.w*AirToFuelFlow)
        
        # set composition of combustion products
        self.out_station.comp = self.current_fluid.burn_comp(self.in_station.comp, cl.Composition(), self.fgrc)
        
        # calc enthalpy of air at entry
        self.in_station.h = self.current_fluid.tp2h(self.in_station.t, self.in_station.p, self.in_station.comp)
        
        # calc enthalpy of combustion products
        self.out_station.h = self.current_fuel.h_fuel(self.wf, self.eta_b, self.in_station.w, self.in_station.h, self.out_station.comp)
        
        # calculate outlet temperature
        self.out_station.t = self.current_fluid.hp2t(self.out_station.h, self.in_station.p, self.out_station.comp)
        
        # calculate outlet mass flow
        self.out_station.w = self.in_station.w + self.wf/AirToFuelFlow*self.eta_b
    
    def reset(self):
        self=Combustor()