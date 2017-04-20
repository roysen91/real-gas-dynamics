'''
Created on 06.02.2017

@author: roysonntag
'''
import classes as cl
from constants import *


class Fuel:

    def __init__(self,hv=None,tf=None):
        self.hv = 10337.74*4184  # [CHU/lbm] to [J/kg]fuel heating value
        self.tf = 382.5948
        self.t_ref = 298.15 # reference temperature for fuel
        self.burn_model = 'CEA'
    
    def h_fuel(self,wf,eta_b,w0,h0,comp):
        if self.burn_model == 'CEA':
        
            fuel_comp = cl.Composition(comp_def='Jet-A1')
            my_fluid = cl.CeaFluid()
            
            enthalpy_fuel = my_fluid.tp2h(self.tf, StandartPressure, fuel_comp)
            
            # fuel to air mass ratio
            fmr = wf*eta_b/(w0*AirToFuelFlow)
            
            # enthalpy combustion products
            return (h0+enthalpy_fuel*fmr)/(1+fmr)
        
        elif self.burn_model == 'RRD':
            comp_air = cl.Composition('Air')
            comp_stoi = cl.Composition('Air')
            comp_stoi.FAR = FARStoi
            my_fluid = cl.RRDFluid()
            
            h_air = my_fluid.tp2h(self.t_ref, StandartPressure, comp_air)
            h_stoi = my_fluid.tp2h(self.t_ref, StandartPressure, comp_stoi)
            
            heating_corr = h_stoi+(h_stoi-h_air)/comp_stoi.FAR
            
            # fuel to air mass ratio
            fmr = wf*eta_b/(w0*AirToFuelFlow)
            # corrected fuel heating value
            enthalpy_fuel = self.hv + heating_corr # +CP*(self.tf-self.t_ref)
            # enthalpy combustion products
            return (h0+enthalpy_fuel*fmr)/(1+fmr)
            

        