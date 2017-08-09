'''
Created on 17.11.2016

@author: roysonntag
'''
from math import *
from constants import *
import numpy as np

class EquationOfState:
    def __init__(self):
        self.a=0
        self.b=0
        self.f=0
        self.alpha=0
    def get_molare_volume(self,t,p,comp):
        min_val = MinMolVol
        iteration = 0
        norm = 1
        vol_mol_i = StartValueMolareVolume
        while iteration<MaxIterations and norm > MinDeviation:
            # check if molare volume is in range
            if vol_mol_i < min_val: vol_mol_i = min_val
            
            f1 = (self.get_pressure(t,vol_mol_i,comp)-p)
            jacobi1 = self.get_derivative_vol_mol(t,vol_mol_i,comp)
            
            # get vol_mol step
            delta_vol_mol = -f1/jacobi1
            vol_mol_i +=delta_vol_mol
            #norm of first dimension
            norm = sqrt(delta_vol_mol**2)            
            iteration +=1
        if iteration==40:
            print('Failure to converge in get_molare_volume for t=',str(t),' and p=',str(p))
        return vol_mol_i
        
            
            
    #def set_coeffs(self,comp):
        # calc real gas constants    
       # elif self.eos_model is 'BER':
        #    self.a = (9/8)*comp.V_mol_c*UnivGasConstant*comp.t_crit**2
         #   self.b = comp.V_mol_c/3
        #elif self.eos_model is 'CLA':
        #    self.a = (27/20)*UnivGasConstant*comp.t_crit**2*comp.V_mol_c
         #   self.b = comp.V_mol_c/5
        #else:
         #   self.a=0
          #  self.b=0

class VanDerWaals(EquationOfState):
    def __init__(self):
        EquationOfState.__init__(self)
    def set_coeffs(self, comp):
        self.a = (27/64)*UnivGasConstant**2*comp.t_crit**2/comp.p_crit
        self.b = (1/8)*UnivGasConstant*comp.t_crit/comp.p_crit
    def get_pressure(self,t,v_mol,comp):
        return UnivGasConstant*t/(v_mol-self.b)-self.a/(v_mol*v_mol)
    def get_derivative_vol_mol(self,t,v_mol,comp):
        return -UnivGasConstant*t/(v_mol-self.b)**2+2*self.a/v_mol**3
    def get_derivative_t(self,t,v_mol):
        return UnivGasConstant/(v_mol-self.b)
    def get_deviation(self,t,p,comp,prop):
        if prop == 'cp':
            v_mol = self.get_molare_volume(t, p,comp)
            return -t/comp.mol_wgt*self.get_derivative_t(t, v_mol)**2*self.get_derivative_vol_mol(t, v_mol, comp)**-1-UnivGasConstant/comp.mol_wgt

class Soave(EquationOfState):
    def __init__(self):
        EquationOfState.__init__(self)
    def set_coeffs(self, comp):
        self.a = (0.42747*UnivGasConstant**2*comp.t_crit**2/comp.p_crit)
        self.b = 0.08664*UnivGasConstant*comp.t_crit/comp.p_crit
        self.f = 0.480 + 1.574 *comp.omega - 0.176 * comp.omega**2
    def get_alpha(self,t,comp):
        return 1+self.f*(1-sqrt(t/comp.t_crit))**2
    def get_pressure(self,t,v_mol,comp):
        return UnivGasConstant*t/(v_mol-self.b)-(self.a*self.get_alpha(t, comp))/(v_mol*(v_mol+self.b))
    def get_derivative_vol_mol(self,t,v_mol,comp):
        return -UnivGasConstant*t/(v_mol-self.b)**2+(self.a*self.get_alpha(t, comp)*(self.b+2*v_mol))/(v_mol**2*(self.b+v_mol)**2)
    def get_derivative_t(self,t,v_mol):
        return UnivGasConstant/(v_mol-self.b)
    def get_deviation(self,t,p,comp,prop):
        if prop == 'cp':
            v_mol = self.get_molare_volume(t, p, comp)
            cp_dev = -t/comp.mol_wgt*self.get_derivative_t(t, v_mol)**2*self.get_derivative_vol_mol(t, v_mol, comp)**-1-UnivGasConstant/comp.mol_wgt
            return cp_dev
        elif prop == 'h':
            v_mol = self.get_molare_volume(t, p, comp)
            h_dev = p*v_mol/comp.mol_wgt-UnivGasConstant*t/comp.mol_wgt+(self.get_alpha(t, comp)*self.a/(self.b*comp.mol_wgt))*log(v_mol/(v_mol-self.b))
            return h_dev
        elif prop == 's':
            v_mol = self.get_molare_volume(t, p, comp)
            s_dev = (UnivGasConstant/comp.mol_wgt)*log((p*(v_mol-self.b))/(UnivGasConstant*t))
            return s_dev

class Species:
    t_crit=0
    p_crit=0
    vol_mol_crit=0
    mol_wgt=0
    omega=0
    name=''
    
    def __init__(self,name):
        self.cea_const=[]
        self.bucker_const=[]
        self.name = name
        self.r_spec = 0

    
        species_file = open('data/species.txt','r')
        for line in species_file:
            line=line.split()
            if line[0] == self.name:
                self.t_crit         = float(line[1])
                self.p_crit         = float(line[2])
                self.vol_mol_crit   = float(line[3])
                self.mol_wgt        = float(line[4])/1000
                self.omega          = float(line[5])
                break
        # will only exec if for loop finishes (end of file hit)
        else:
            print('Could not find Species:',self.name,' in species.txt!')
        # set specific gas constant
        self.r_spec = UnivGasConstant/self.mol_wgt
        # load CEA coeffs from file
        cea_file = open('data/cea_constants.txt','r')
        for line in cea_file:
            line=line.split()
            if line[0] == self.name:
                for value in range(1,len(line)):
                    self.cea_const.append(float(line[value]))
                break
        # will only exec if for loop finishes (end of file hit)
        else:
            print('Could not find CEA coeffs for:',self.name,' in cea_constants.txt!')
        # load Bucker coeffs from file
        bucker_file = open('data/bucker_constants.txt','r')
        for line in bucker_file:
            line=line.split(',')
            if line[0] == self.name:
                for value in range(1,len(line)):
                    self.bucker_const.append(float(line[value]))
                break
        # will only exec if for loop finishes (end of file hit)
        else:
            # there are no bucker constants for jet-a
            pass
            #print('Could not find Bucker coeffs for:',self.name,' in bucker_constants.txt!')

class Composition(Species):

    def __init__(self,name='new_comp',comp_def=None):
        self.structure=dict()       # dictionary holding vol fraction and name of each species
        self.mass_struct = dict()   # dict holding each species mass fraction
        self.rel_wgt=dict()       # dictionary holding mass fraction and name of each species
        self.r_spec=0    # specific gas constant  
        self.FAR = 0
        self.WGR_mass = 0    # Water to dry Gas mass ratio
        self.WGR_mole = 0    # Water to dry Gas mole ratio
        # set composition to standard dry air
        if comp_def == None or comp_def == 'Air':
            self.name = 'Air'
            self.structure = {  Species('O2'): 0.209476, 
                                Species('CO2'): 0.000319,
                                Species('AR'): 0.0093650,
                                Species('N2'): 0.78084}
        # comp is Jet-A1 fuel
        elif comp_def == 'Jet-A1':
            self.name = 'Jet-A1(L)'
            self.FAR = 1
            self.structure = {  Species(self.name): 1} 
        # comp is custom
        else:
            self.name = name
            self.structure=comp_def
        
        self.set_gas_properties()
        
    def set_eos(self,eos_name):
        self.eos = EquationOfState(self,eos_name)
        
    def set_gas_properties(self):
        self.t_crit         = 0
        self.p_crit         = 0
        self.vol_mol_crit   = 0
        self.mol_wgt        = 0
        self.omega          = 0
        
        for sp,fraction in self.structure.items():
            self.t_crit         += sp.t_crit*fraction
            self.p_crit         += sp.p_crit*fraction
            self.vol_mol_crit   += sp.vol_mol_crit*fraction
            self.mol_wgt        += sp.mol_wgt*fraction
            self.omega          += sp.omega*fraction
        #calc relative weight of each species
        for sp,rel_vol in self.structure.items():
            self.mass_struct[sp.name] = rel_vol*(sp.mol_wgt/self.mol_wgt)
        self.r_spec = UnivGasConstant/self.mol_wgt
    
    # returns list of comp fractions
    def get_fraction(self,req_fraction):
        for sp,fraction in self.structure.items():
            if sp.name == req_fraction:
                return fraction
        return 0
    def get_density(self,t,p):
        return (UnivGasConstant*t)/(p*self.mol_wgt)
            
    def set_WGR(self,WGR_mass):
        if WGR_mass == 0:
            #do nothing
            return
        self.WGR_mass = WGR_mass
        water_obj = Species('H2O')
        dry_air_obj = Composition('Air')

        # number of moles of water
        # n_dry_air = 1
        n_water = 1/(1+(water_obj.mol_wgt/dry_air_obj.mol_wgt)*(1/WGR_mass))
        # set molare Water to Gas Ratio
        self.WGR_mole = WGR_mass*(dry_air_obj.mol_wgt/water_obj.mol_wgt)
        # new mole fraction of compositiion
        total = 0
        # set fraction of water
        self.structure[water_obj] = n_water/(1+n_water)
        for sp,frac in self.structure.items():
            new_frac = (frac)/(1+n_water)
            self.structure[sp] = new_frac
            total += new_frac
        if abs(1-total)>0.01:
            print('In set_WGR total of species fractions does not sum up to 1')
        else:
            #set new gas properties
            self.set_gas_properties()

class Fluid():
    # maximal number of iteration for konvergance
    max_iter = 40
    # termination criterion for newton method
    min_dev = 1e-7
    def __init__(self,eos_model='VDW'):
        if eos_model == 'VDW':
            self.eos = VanDerWaals()
        elif eos_model == 'SOA':
            self.eos = Soave()
    def set_eos(self,eos_name):
        self.eos = EquationOfState(self,eos_name)   
    def tp2density(self,t,p,comp):
        comp.get_density(t,p)
    def tp2h(self,t,p,comp,mode='ideal'):
        if mode == 'ideal':
            return self.ideal(t,p,comp,'h')
        elif mode == 'real':
            return self.real(t,p,comp,'h')
    def tp2s(self,t,p,comp,mode='ideal'):
        if mode == 'ideal':
            return self.ideal(t,p,comp,'s')
        elif mode == 'real':
            return self.real(t,p,comp,'s')
    def tp2cp(self,t,p,comp,mode='ideal'):
        if mode == 'ideal':
            return self.ideal(t,p,comp,'cp')
        elif mode == 'real':
            return self.real(t,p,comp,'cp')
    def real(self,t,p,comp,prop):
        return self.ideal(t, p, comp, prop)+self.eos.get_deviation(t,p,comp,prop)
    def tp2kappa(self,t,p,comp,mode='ideal'):
        cp = self.tp2cp(t,p,comp,mode)
        return cp/(cp-comp.r_spec)
    def tp2kappa_eos(self,t,p,comp,mode='ideal'):
        v_mol = self.eos.get_molare_volume(t, p, comp)
        # get expansivity
        alpha_v = - self.eos.get_derivative_t(t, v_mol)/(v_mol*self.eos.get_derivative_vol_mol(t, v_mol, comp))
        # get isothermal compressibility
        K_isotherm = -1/(v_mol*self.eos.get_derivative_vol_mol(t, v_mol, comp))
        # get spec. heat in J/mol*K
        cp = self.tp2cp(t, p, comp, mode)*comp.mol_wgt
        # get isentrop compressibility
        K_isentrop = K_isotherm-(v_mol*t*alpha_v**2)/cp
        # calc kappa
        kappa = K_isotherm/K_isentrop
        
        return kappa
        
    def tp2poly_exponent(self,eta_poly,t,p,comp):
        kappa = self.tp2kappa_ideal(t,p,comp);
        return 1/(1-((kappa-1)/kappa)*(1/eta_poly))
    def isentropic_compression(self,t_in,p_in,p_out,comp,mode='ideal',eta_isentropic=1):
        # kappa at inlet conditions
        kappa_in = self.tp2kappa(t_in,p_in,comp,mode)
        # first guess of tout with constant kappa
        t_out = t_in *(((p_out/p_in)**((kappa_in-1)/kappa_in)-1)*(1/eta_isentropic)+1)
                
        i=0
        norm=1
        dt=0
        while i < self.max_iter and norm  > self.min_dev:
            # kappa at outlet conditions
            kappa_out = self.tp2kappa(t_out+dt,p_out,comp,mode)
            # use arathmetic mean for kappa
            kappa = (kappa_in+kappa_out)/2
            # calc new t_out with new kappa
            t_out_n = t_in *((((p_out/p_in)**((kappa-1)/kappa)-1)/eta_isentropic)+1)
            # new delta in temperature
            dt = t_out_n-t_out
            # t_out for next iteration
            t_out += dt
            # norm in first dimension
            norm = sqrt(dt*dt)
            
            i+=1
        if i==40:
            print('Failure to converge in isentropic_compression for t_in=',str(t_in),' and PR=',str(p_out/p_in))  
        return t_out
    def polytropic_compression(self,eta_poly,t_in,p_in,p_out,comp):
        poly_exp = self.tp2poly_exponent(eta_poly,t_in,p_in,comp)
        t_out = t_in *(p_out/p_in)**((poly_exp-1)/poly_exp)
        return t_out
    def burn_comp(self,gas,fuel,fmr):
        pass
    def hp2t(self,h,p,comp):
        t_n = 1000
        i = 0
        norm = 1
        dt = 0
        
        # use newton rephson method 
        while i < self.max_iter and norm  > self.min_dev:
            if t_n < MinTemp: t_n = MinTemp
            
            # declare next step value
            t0 = t_n *1.001
            
            f1 = (self.tp2h(t_n, p, comp)-h)
            jakobi = (self.tp2h(t_n, p, comp)-self.tp2h(t0, p, comp))/(t_n-t0)
            
            # new delta in temperature
            dt = -f1/jakobi
            
            # new start value
            t_n += dt
            # norm in first dimension
            norm = sqrt(dt*dt)
            
            i+=1
        if i==40:
            print('Failure to converge in HP2T for h=',str(h),' and p=',str(p))
        return t_n
    def sp2t(self,s,p,comp):
        t_n = 1000
        i = 0
        norm = 1
        dt = 0
        
        # use newton rephson method 
        while i < self.max_iter and norm  > self.min_dev:
            if t_n < MinTemp: t_n = MinTemp
            
            # declare next step value
            t0 = t_n *1.001
            
            f1 = (self.tp2s(t_n, p, comp)-s)
            jakobi = (self.tp2s(t_n, p, comp)-self.tp2s(t0, p, comp))/(t_n-t0)
            
            # new delta in temperature
            dt = -f1/jakobi
            
            # new start value
            t_n += dt
            # norm in first dimension
            norm = sqrt(dt*dt)
            
            i+=1
        if i==40:
            print('Failure to converge in SP2T for s=',str(s),' and p=',str(p))
        return t_n
    def get_p_sat(self,t):
        '''calculates saturation pressure for air at given temperature
        INPUT:     float t        temperature           [K]
        OUTPUT:    float p_sat    saturation pressure    [Pa]
        '''
        t_cel = t-273.15
        if (t_cel<=80):
            p_sat_1 = 133.32185 * exp(18.46882-(3827.77+218140.0/t)/t)
        elif t_cel >=70:
            M1 = -7.691234564
            M2 = -26.08023696
            M3 = -168.1706546
            M4 = 64.23285504
            M5 = -118.9646225
            M6 = 4.167117320
            M7 = 20.97506760
            M8 = 1000000000.0
            M9 = 6.0
    
            theta = t / 647.3
    
            m_sat_b = M1 * (1.0 - theta)+ M2 * pow(1.0 - theta, 2)+ M3 * pow(1.0 - theta, 3)+ M4 * pow(1.0 - theta, 4)+ M5 * pow(1.0 - theta, 5)
    
            # beta = saturation beta at theta
            beta = exp((((1.0 / theta) * m_sat_b) / (1.0 + (M6 * (1.0 - theta))+(M7 * (pow(1.0 - theta, 2)))))
                    - ((1.0 - theta) / ((M8 * (pow(1.0 - theta, 2))) + M9)))
    
            p_sat_2 = beta * 22120000.0
        # use Psat1 or Psat2 or combination depending on t_cel
        if (t_cel <= 70.): 
            return p_sat_1
        elif (t_cel > 70. and t_cel < 80.): 
            return p_sat_1 + (t_cel - 70) / 10 * (p_sat_2 - p_sat_1);
        else: 
            return p_sat_2;
    def set_humidity(self,t,p,comp,WGR_mass=None,RH=None,SH=None):
        ''' sets WGR_mass of comp depending in humidity input mode and value
            INPUT:  float    t        static temperature
                    float    p        static pressure
                    obj      comp     Composition object
                    float    WGR_mass water to dry air mass ratio
                    float    RH       relative humidity
                    float    SH       specific humidity
            OUTPUT:
                    void    
        '''
        p_sat = self.get_p_sat(t)
        t_dat = 243.15
        air_obj = Composition('Air')
        water_obj = Species('H2O')
        # avoid division by zero
        if p-p_sat==0:        
            print('Saturation pressure equals static pressure. Humidity calculation not possible due to division by zero!')
        else:  
            WGR_max = (air_obj.r_spec/water_obj.r_spec)*p_sat/(p-p_sat)       
            if WGR_mass!=None and RH==None and SH==None:
                # set WGR directly
                comp.set_WGR(WGR_mass)
            elif WGR_mass==None and RH!=None and SH==None:
                if p_sat < p and t >= t_dat:
                    if RH <= 0:
                        comp.set_WGR(0)
                    else:
                        if RH <= 1:
                            comp.set_WGR((air_obj.r_spec/water_obj.r_spec)*p_sat*RH/(p-p_sat*RH))
                        else:
                            print('Relative humidity invalid: RH = '+str(RH))
                        if comp.WGR_mass>WGR_max:
                            print('Water to gas ratio: WGR = '+str(WGR_mass)+'> Max value:',+WGR_max,'reset to relative humidity = 1.0')
                            comp.set_WGR(WGR_max)
                else:
                    # WGR limited by saturation
                    if RH>0:
                        if t<t_dat:
                            print('Ambient Temperature: T='+str(t)+'<'+str(t_dat)+', humidity set to zero because equation for p_sat not valid')  
                            comp.set_WGR(0.0) 
                        else:
                            comp.set_WGR(1e30) 
                    else:  
                        comp.set_WGR(0.0) 
            elif WGR_mass==None and RH==None and SH!=None:
                print('No method for specific humidity in place')    

class RRDFluid(Fluid):
    def __init__(self,eos_model='VDW'):
        Fluid.__init__(self,eos_model)
    def set_eos(self,eos_name):
        self.eos = EquationOfState(self,eos_name)   
    def ideal(self,t,p,comp,prop):
        if prop is 'h':
            x = 15.682994156168 * comp.FAR / (1. + comp.FAR)
            # low temperature range
            if (t <= 1000.0):  
                # enthalpy for air and products of combustion
                h_ges = (256.800991 + 81.3598929 * x + 
                        t * (5.978055E-1 + 4.015965E-1 * x +
                        t * (9.78422E-5 - 1.947137E-3 * x +
                        t * (-3.950606E-7 + 5.391760E-6 * x +
                        t * (7.991019E-10 - 8.424275E-9 * x +
                        t * (-7.154873E-13 + 7.595285E-12 * x +
                        t * (2.973926E-16 - 3.683035E-15 * x +
                        t * (-4.544470E-20 + 7.441810E-19 * x))))))))
                # adjust enthalpy for extra water content
                if (comp.WGR_mass != 0.0):
                    h_water = (1738.925309 +
                            t * (1.165116 +
                            t * (-2.357109E-4 +
                            t * (5.224504E-7 +
                            t * (-3.136553E-10 +
                            t * (7.697308E-14))))))
                    h_ges = ((comp.WGR_mass * h_water) + h_ges) / (1.0 + comp.WGR_mass)
            # high temperature range
            elif (t > 1000.0):
                # enthalpy for air and products of combustion
                h_ges = (269.6606 + 115.2169 * x +
                        t * (5.302190E-1 + 7.855338E-3 * x +
                        t * (1.170436E-4 + 4.526546E-5 * x +
                        t * (-2.649052E-8 - 9.34204E-9 * x +
                        t * (2.517911E-12 + 7.397349E-13 * x)))))
                # adjust enthalpy for extra water content
                if (comp.WGR_mass != 0.0):
                    h_water = (1847.551754 +
                            t * (7.662467E-1 +
                            t * (4.016958E-4 +
                            t * (-6.583356E-8 +
                            t * (4.437895E-12)))))
                    h_ges = ((comp.WGR_mass * h_water) + h_ges) / (1.0 + comp.WGR_mass)
            return h_ges * RRDSpecHeatConversion * SpecHeatConversion

        elif prop is 's':
            x = 15.682994156168 * comp.FAR / (1. + comp.FAR)
            if (t <= 1000.0):  # low temperature range
                # entropy function for air and combustion products
                s_ges = (0.640479618837074 - 1.61852820336 * x + 
                        (5.978055E-1 + 4.015965E-1 * x) * log(t) +
                        t * (1.956844E-4 - 3.894273E-3 * x +
                        t * (-5.925909E-7 + 8.087640E-6 * x +
                        t * (1.065469E-9 - 1.123237E-8 * x +
                        t * (-8.943591E-13 + 9.49411E-12 * x +
                        t * (3.568711E-16 - 4.419641E-15 * x +
                        t * (-5.301881E-20 + 8.682111E-19 * x)))))))
                # adjust entropy function for extra water content
                if (comp.WGR_mass != 0.0):
                    s_water = (-0.1853961196989 +
                            (1.165116 * log(t) +
                            t * (-4.714218E-4 +
                            t * (7.836755E-7 +
                            t * (-4.182070E-10 +
                            t * (9.621635E-14))))))
                    s_ges = ((comp.WGR_mass * s_water) + s_ges) / (1.0 + comp.WGR_mass)
            elif (t > 1000.0):  # high temperature range
                # entropy function for air and combustion products /
                s_ges = (0.987720 - 0.072510 * x +
                        (5.302190E-1 + 7.855338E-3 * x) * log(t) +
                        t * (2.340672E-4 + 9.05509E-5 * x +
                        t * (-3.973577E-8 - 1.401306E-8 * x +
                        t * (3.354881E-12 + 9.88647E-13 * x))))
                # adjust entropy function for extra water content /
                if (comp.WGR_mass != 0.0):  # water
                    s_water = (1.84960 +
                            (7.662467E-1 * log(t) +
                            t * (8.033916E-4 +
                            t * (-9.875035E-8 +
                            t * (5.917193E-12)))))
                    s_ges = ((comp.WGR_mass * s_water) + s_ges) / (1.0 + comp.WGR_mass)
            # Specific Entropy
            return s_ges * RRDSpecHeatConversion* SpecHeatConversion  - (UnivGasConstant/comp.mol_wgt) * (log(p) - log(StandartPressure))

        elif prop is 'cp':
            x = 15.682994156168 * comp.FAR / (1. + comp.FAR)

            # specific heat for air and products of combustion
            if (t <= 1000.0): # low temperature range
                cp_ges = (5.978055E-1 + 4.015965E-1 * x +
                        t * (2.0 * 9.78422E-5 - 2.0 * 1.947137E-3 * x +
                        t * (-3.0 * 3.950606E-7 + 3.0 * 5.391760E-6 * x +
                        t * (4.0 * 7.991019E-10 - 4.0 * 8.424275E-9 * x +
                        t * (-5.0 * 7.154873E-13 + 5.0 * 7.595285E-12 * x +
                        t * (6.0 * 2.973926E-16 - 6.0 * 3.683035E-15 * x +
                        t * (-7.0 * 4.544470E-20 + 7.0 * 7.441810E-19 * x)))))))
                # adjust specific heat for extra water content
                if (comp.WGR_mass != 0.0): # water
                    cp_water = (1.165116 +
                            t * (-2.0 * 2.357109E-4 +
                            t * (3.0 * 5.224504E-7 +
                            t * (-4.0 * 3.136553E-10 +
                            t * (5.0 * 7.697308E-14)))))
                    cp_ges = (comp.WGR_mass * cp_water + cp_ges) / (1.0 + comp.WGR_mass)
            # specific heat for air and products of combustion
            elif  (t > 1000):  # high temperature range
                cp_ges = (5.302190E-1 + 7.855338E-3 * x +
                        t * (2.0 * 1.170436E-4 + 2.0 * 4.526546E-5 * x +
                        t * (-3.0 * 2.649052E-8 - 3.0 * 9.34204E-9 * x +
                        t * (4.0 * 2.517911E-12 + 4.0 * 7.397349E-13 * x))))
                # adjust specific heat for extra water content
                if (comp.WGR_mass != 0.0): # water
                    cp_water = (7.662467E-1 + 
                            t * (2.0 * 4.016958E-4 +
                            t * (-3.0 * 6.583356E-8 +
                            t * (4.0 * 4.437895E-12))))
                        
                    cp_ges = (comp.WGR_mass * cp_water + cp_ges) / (1.0 + comp.WGR_mass)
            return cp_ges * RRDSpecHeatConversion * SpecHeatConversion
    def burn_comp(self, gas, fuel, fmr):
        Fluid.burn_comp(self, gas, fuel, fmr)
        
        comp = Composition('burn_prod')
        comp.FAR = gas.FAR + fmr *(1+gas.FAR)*(1+gas.WGR_mass)
        comp.WGR_mass = gas.WGR_mass/(1+fmr*(1+gas.WGR_mass))
        
        return comp

class BuckerFluid(Fluid):
    def __init__(self,eos_model='VDW'):
        Fluid.__init__(self,eos_model)
        self.t0 = 273.15
        # power for later use
        self.b = [0.00, -1.50, -1.25, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1.0]
    def ideal(self,t,p,comp,prop):
        if prop is 'h':
            h_ges = 0
            h_species = 0
            for species,fraction in comp.structure.items():
                for i in range(10):
                    h_species += species.bucker_const[i+2]*self.t0/(self.b[i]+1)*(t/self.t0)**(self.b[i]+1)
                h_ges += fraction*(h_species+species.bucker_const[0])/comp.mol_wgt
                h_species = 0
            return h_ges

        elif prop is 's':
            s_ges = 0
            s_species = 0
            s_mix_k = 0
            for species,fraction in comp.structure.items():
                for i in range(3,12):
                    s_species += species.bucker_const[i]/self.b[i-2]*(t/self.t0)**self.b[i-2]
                s_mix_k += (fraction*(species.bucker_const[1]-UnivGasConstant*log(p/StandartPressure)+species.bucker_const[2]*log(t/self.t0)+s_species)
                           -(UnivGasConstant/comp.mol_wgt)*fraction*log(fraction))
                s_species = 0
            return s_mix_k -(UnivGasConstant/comp.mol_wgt)*log(p/StandartPressure)

        elif prop is 'cp':
            cp_ges = 0
            cp_species = 0
            for species,fraction in comp.structure.items():
                for i in range(10):
                    cp_species += species.bucker_const[i+2]*(t/self.t0)**self.b[i]
                cp_ges += fraction*(cp_species/comp.mol_wgt)
                cp_species = 0
            return cp_ges

class CeaFluid(Fluid):
    def __init__(self,eos_model='VDW'):
        Fluid.__init__(self,eos_model)
    def ideal(self,t,p,comp,prop):
        if prop is 'h':
            h_ges=0
            if t<=1000:
                # sum up enthalpy of every single species
                for sp,fraction in comp.structure.items():
                    # sum up using species specific CEA constants                    
                    h_ges+=(fraction*UnivGasConstant*t/(comp.mol_wgt))*(-sp.cea_const[0]/t**2                                                                            
                                                                      +sp.cea_const[1]*log(t)/t
                                                                      +sp.cea_const[2]
                                                                      +sp.cea_const[3]*t/2
                                                                      +sp.cea_const[4]*t**2/3
                                                                      +sp.cea_const[5]*t**3/4
                                                                      +sp.cea_const[6]*t**4/5
                                                                      +sp.cea_const[7]/t)
            else:
                # sum zp enthalpy of every single species
                for sp,fraction in comp.structure.items():
                    # sum up using species specific CEA constants     
                    h_ges+=(fraction*UnivGasConstant*t/(comp.mol_wgt))*(-sp.cea_const[9]/t**2   
                                                                      +sp.cea_const[10]*log(t)/t
                                                                      +sp.cea_const[11]
                                                                      +sp.cea_const[12]*t/2
                                                                      +sp.cea_const[13]*t**2/3
                                                                      +sp.cea_const[14]*t**3/4
                                                                      +sp.cea_const[15]*t**4/5
                                                                      +sp.cea_const[16]/t)
            return h_ges
        elif prop is 's':
            s_ges=0
            s_mix=0
            if t<=1000:
                # sum up enthalpy of every single species
                for sp,fraction in comp.structure.items():
                    # sum up using species specific CEA constants
                    s_ges+=(fraction*UnivGasConstant*sp.mol_wgt/(comp.mol_wgt*comp.mol_wgt))*(-sp.cea_const[0]/(2*t**2)
                                                                      -sp.cea_const[1]/t
                                                                      +sp.cea_const[2]*log(t)
                                                                      +sp.cea_const[3]*t
                                                                      +sp.cea_const[4]*t**2/2
                                                                      +sp.cea_const[5]*t**3/3
                                                                      +sp.cea_const[6]*t**4/4
                                                                      +sp.cea_const[8])
                    s_mix -= fraction*(UnivGasConstant/comp.mol_wgt)*log(fraction)
                s_mix -= (UnivGasConstant/comp.mol_wgt)*(log(p)-log(StandartPressure))
            else:
                # sum up enthalpy of every single species
                for sp,fraction in comp.structure.items():
                    # sum up using species specific CEA constants
                    s_ges+=(fraction*UnivGasConstant*sp.mol_wgt/(comp.mol_wgt*comp.mol_wgt))*(-sp.cea_const[9]/(2*t**2)
                                                                      -sp.cea_const[10]/t
                                                                      +sp.cea_const[11]*log(t)
                                                                      +sp.cea_const[12]*t
                                                                      +sp.cea_const[13]*t**2/2
                                                                      +sp.cea_const[14]*t**3/3
                                                                      +sp.cea_const[15]*t**4/4
                                                                      +sp.cea_const[17])
                    s_mix -= fraction*(UnivGasConstant/comp.mol_wgt)*log(fraction)
                s_mix -= (UnivGasConstant/comp.mol_wgt)*(log(p)-log(StandartPressure))
            return s_ges+s_mix
        elif prop is 'cp':
                cp_ges=0
                if t<=1000:
                    # sum up enthalpy of every single species
                    for sp,fraction in comp.structure.items():
                        # sum up using species specific CEA constants
                        cp_ges+=(fraction*UnivGasConstant/(comp.mol_wgt))*(+sp.cea_const[0]/t**2
                                                                          +sp.cea_const[1]/t
                                                                          +sp.cea_const[2]
                                                                          +sp.cea_const[3]*t
                                                                          +sp.cea_const[4]*t**2
                                                                          +sp.cea_const[5]*t**3
                                                                          +sp.cea_const[6]*t**4)
                else:
                    # sum up enthalpy of every single species
                    for sp,fraction in comp.structure.items():
                        # sum up using species specific CEA constants
                        cp_ges+=(fraction*UnivGasConstant/(comp.mol_wgt))*(+sp.cea_const[9]/t**2
                                                                          +sp.cea_const[10]/t
                                                                          +sp.cea_const[11]
                                                                          +sp.cea_const[12]*t
                                                                          +sp.cea_const[13]*t**2
                                                                          +sp.cea_const[14]*t**3
                                                                          +sp.cea_const[15]*t**4)
                return cp_ges
    def burn_comp(self, gas, fuel, fmr):
        Fluid.burn_comp(self, gas, fuel, fmr)
        
        # coeffs for Air composition
        c1 = gas.get_fraction('N2')/gas.get_fraction('O2')
        c2 = gas.get_fraction('CO2')/gas.get_fraction('O2')
        c3 = gas.get_fraction('AR')/gas.get_fraction('O2')
        
        dummy = Composition()
        b1 = dummy.get_fraction('N2')/dummy.get_fraction('O2')
        b2 = dummy.get_fraction('CO2')/dummy.get_fraction('O2')
        b3 = dummy.get_fraction('AR')/dummy.get_fraction('O2')
        
        # Jet A C12H23 used here
        mol_of_C = 12
        mol_of_H = 23
        
        eq_ratio = fmr/FARStoi
        fuel_param = mol_of_C+mol_of_H/4
        
        mol_of_water = gas.WGR_mass *((fuel_param*(Species('O2').mol_wgt+c1*Species('N2').mol_wgt+c2*Species('CO2').mol_wgt+c3*Species('AR').mol_wgt))/(Species('H2O').mol_wgt))
        
        coeffs = np.zeros(5)
        
        coeffs[0]=eq_ratio*mol_of_H/2+mol_of_water  # H2O
        coeffs[1]=eq_ratio*mol_of_C+fuel_param*c2   # CO2
        coeffs[2]=fuel_param*(1-eq_ratio)           # O2
        coeffs[3]=fuel_param*c1                     # N2
        coeffs[4]=fuel_param*c3                     # AR
        
        # total number of moles per species
        nt = np.sum(coeffs)
        
        # deviation from value one
        dev = abs(np.sum(coeffs/nt)-1)
        
        if dev > 0.01:
            print('In burn_comp sum of species fractions don not sum up to 1')
        else:
            combustion_prod = Composition(comp_def={Species('H2O'): coeffs[0]/nt,Species('CO2'): coeffs[1]/nt, Species('O2'): coeffs[2]/nt, Species('N2'): coeffs[3]/nt, Species('AR'): coeffs[4]/nt})
        
        return combustion_prod