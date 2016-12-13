'''
Created on 17.11.2016

@author: roysonntag
'''
from math import *
UnivGasConstant= 8.3143 # [J/mol*K]
StandartPressure = 101325 # [Pa]

class EOS:
    
    a=0
    b=0
    f=0
    alpha=0
    eos_model='VDW'
    #UnivGasConst=8.3143       # [J/mol*K] universal gas constant
       
    def __init__(self,comp=None,eos_model=None):
        if comp is None:
            comp=Composition('Air')
        #   overwirte default eos 
        self.eos_model = eos_model
        # calc real gas constants    
        if self.eos_model is 'VDW':
            self.a = (27/64)*UnivGasConstant**2*comp.t_crit**2/comp.p_crit;
            self.b = (1/8)*UnivGasConstant*comp.t_crit/comp.p_crit;
        elif self.eos_model is 'BER':
            self.a = (9/8)*comp.V_mol_c*UnivGasConstant*comp.t_crit**2;
            self.b = comp.V_mol_c/3;
        elif self.eos_model is 'CLA':
            self.a = (27/20)*UnivGasConstant*comp.t_crit**2*comp.V_mol_c;
            self.b = comp.V_mol_c/5;
        elif self.eos_model is 'SOA':
            self.a = (0.42747*UnivGasConstant**2*comp.t_crit**2/comp.p_crit);
            self.b = 0.08664*UnivGasConstant*comp.t_crit/comp.p_crit;
            self.f = 0.480 + 1.574 *comp.omega - 0.176 * comp.omega**2;
        else:
            self.a=0
            self.b=0
        
    
class Species:
    name=''
    t_crit=0
    p_crit=0
    vol_mol_crit=0
    mol_wgt=0
    omega=0
  
    
    def __init__(self,name):
        self.cea_const=[]
        self.bucker_const=[]
        self.name = name
        species_file = open('species.txt','r')
        for line in species_file:
            line=line.split()
            if line[0] == self.name:
                self.t_crit         = float(line[1])
                self.p_crit         = float(line[2])
                self.vol_mol_crit   = float(line[3])
                self.mol_wgt        = float(line[4])
                self.omega          = float(line[5])
                break
        # load CEA coeffs from file
        cea_file = open('cea_constants.txt','r')
        for line in cea_file:
            line=line.split()
            if line[0] == self.name:
                for value in range(1,len(line)):
                    self.cea_const.append(float(line[value]))
                break
        # load Bucker coeffs from file
        bucker_file = open('bucker_constants.txt','r')
        for line in bucker_file:
            line=line.split(',')
            if line[0] == self.name:
                for value in range(1,len(line)):
                    self.bucker_const.append(float(line[value]))
                break

        
class Composition(Species):
    num=0           # number of species
    structure=dict()       # dictionary holding mass fraction and name of each species
    r_spec=0    # specific gas constant  
    name=''
    def __init__(self,comp_def=None):
        # set composition to standard dry air
        if comp_def is None or comp_def is 'Air':
            self.name = 'Air'
            self.structure = {  Species('O2'): 0.209476, 
                                Species('CO2'): 0.000319,
                                Species('AR'): 0.0093650,
                                Species('N2'): 0.78084}
        # comp is Jet-A1 fuel
        elif comp_def is 'Jet-A1':
            self.name = 'Jet-A1'
            self.structure = {  comp_def: 1} 
        # comp is custom
        else:
            self.name = 'Mixture'
            self.structure=comp_def
        for sp,fraction in self.structure.items():
            self.t_crit         += sp.t_crit*fraction
            self.p_crit         += sp.p_crit*fraction
            self.vol_mol_crit   += sp.vol_mol_crit*fraction
            self.mol_wgt        += (sp.mol_wgt/1000)*fraction
            self.omega          += sp.omega*fraction
    def set_eos(self,eos_name):
        self.eos = EOS(self,eos_name)

class Fluid(EOS):
    
    def tp2h(self,t,p,comp):
        return self.ideal(t,p,comp,'h')
        #return self.ideal(t,p,comp,'h')+self.dep(t,p,comp,'h)
    def tp2s(self,t,p,comp):
        return self.ideal(t,p,comp,'s')
    def tp2cp(self,t,p,comp):
        return self.ideal(t,p,comp,'cp')
class BuckerFluid(Fluid):
    def __init__(self):
        self.eos=EOS()
        self.t0 = 273.15
        # power for later use
        self.b = [0.00, -1.50, -1.25, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1.0]
    def set_eos(self,eos_name):
        self.eos = EOS(self,eos_name)   
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
    def __init__(self):
        self.eos=EOS()
    def set_eos(self,eos_name):
        self.eos = EOS(self,eos_name)   
    def ideal(self,t,p,comp,prop):
        if prop is 'h':
            h_ges=0
            if t<=1000:
                # sum up enthalpy of every single species
                for sp,fraction in comp.structure.items():
                    # sum up using species specific CEA constants
                    h_ges+=(fraction*UnivGasConstant*t/comp.mol_wgt)*(-sp.cea_const[0]/t**2
                                                                      +sp.cea_const[1]*log(t)/t
                                                                      +sp.cea_const[2]
                                                                      +sp.cea_const[3]*t/2
                                                                      +sp.cea_const[4]*t**2/3
                                                                      +sp.cea_const[5]*t**3/4
                                                                      +sp.cea_const[6]*t**4/5
                                                                      +sp.cea_const[7]/t)
            else:
                # sum zp enthalpy of every single species
                for name,fraction in comp.structure.items():
                    # create species object
                    sp = Species(name)
                    # sum up using species specific CEA constants
                    h_ges+=(fraction*UnivGasConstant*t/comp.mol_wgt)*(-sp.cea_const[9]/t**2
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
                    s_ges+=(fraction*UnivGasConstant/comp.mol_wgt)*(-sp.cea_const[0]/(2*t**2)
                                                                      -sp.cea_const[1]/t
                                                                      +sp.cea_const[2]*log(t)
                                                                      +sp.cea_const[3]*t
                                                                      +sp.cea_const[4]*t**2/2
                                                                      +sp.cea_const[5]*t**3/3
                                                                      +sp.cea_const[6]*t**4/4
                                                                      +sp.cea_const[8])
                    s_mix+=-fraction*(UnivGasConstant/comp.mol_wgt)*log(fraction)
            else:
                # sum up enthalpy of every single species
                for name,fraction in comp.structure.items():
                    # create species object
                    sp = Species(name)
                    # sum up using species specific CEA constants
                    s_ges+=(fraction*UnivGasConstant/comp.mol_wgt)*(-sp.cea_const[9]/(2*t**2)
                                                                      -sp.cea_const[10]/t
                                                                      +sp.cea_const[11]*log(t)
                                                                      +sp.cea_const[12]*t
                                                                      +sp.cea_const[13]*t**2/2
                                                                      +sp.cea_const[14]*t**3/3
                                                                      +sp.cea_const[15]*t**4/4
                                                                      +sp.cea_const[17])
                s_mix+=-fraction*(UnivGasConstant/comp.mol_wgt)*log(fraction)
            return s_ges+s_mix
        elif prop is 'cp':
                cp_ges=0
                if t<=1000:
                    # sum up enthalpy of every single species
                    for sp,fraction in comp.structure.items():
                        # sum up using species specific CEA constants
                        cp_ges+=(fraction*UnivGasConstant/comp.mol_wgt)*(+sp.cea_const[0]/t**2
                                                                          +sp.cea_const[1]/t
                                                                          +sp.cea_const[2]
                                                                          +sp.cea_const[3]*t
                                                                          +sp.cea_const[4]*t**2
                                                                          +sp.cea_const[5]*t**3
                                                                          +sp.cea_const[6]*t**4)
                else:
                    # sum up enthalpy of every single species
                    for name,fraction in comp.structure.items():
                        # create species object
                        sp = Species(name)
                        # sum up using species specific CEA constants
                        cp_ges+=(fraction*UnivGasConstant/comp.mol_wgt)*(+sp.cea_const[9]/(2*t**2)
                                                                          +sp.cea_const[10]/t
                                                                          +sp.cea_const[11]
                                                                          +sp.cea_const[12]*t
                                                                          +sp.cea_const[13]*t**2
                                                                          +sp.cea_const[14]*t**3
                                                                          +sp.cea_const[15]*t**4)
                return cp_ges
