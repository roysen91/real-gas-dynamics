'''
Created on 17.11.2016

@author: roysonntag
'''
from math import *
UnivGasConstant= 8.3143 # [J/mol*K]

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
        species_file = open('cea_constants.txt','r')
        for line in species_file:
            line=line.split()
            if line[0] == self.name:
                for value in range(1,len(line)):
                    self.cea_const.append(float(line[value]))
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
            self.mol_wgt        += sp.mol_wgt*fraction
            self.omega          += sp.omega*fraction
    def set_eos(self,eos_name):
        self.eos = EOS(self,eos_name)

class Fluid(EOS):
    ideal_method='CEA'  
    def __init__(self):
        self.eos=EOS()
    def set_eos(self,eos_name):
        self.eos = EOS(self,eos_name)   
    def set_ideal_method(self,name):
        # CEA,BUCKER,RRD,BRAUN
        self.ideal_method=name
    def ideal(self,t,p,comp,prop):
        if prop is 'h':
            h_ges=0
            if t<=1000:
                # sum zp enthalpy of every single species
                for sp,fraction in comp.structure.items():
                    # sum up using species specific CEA constants
                    h_ges+=(fraction*UnivGasConstant*t/comp.mol_wgt)*(-sp.cea_const[0]/t**2+sp.cea_const[1]*log(t)/t+sp.cea_const[2]+sp.cea_const[3]*t/2+sp.cea_const[4]*t**2/3+sp.cea_const[5]*t**3/4+sp.cea_const[6]*t**4/5+sp.cea_const[7]/t)
            else:
                # sum zp enthalpy of every single species
                for name,fraction in comp.structure.items():
                    # create species object
                    sp = Species(name)
                    # sum up using species specific CEA constants
                    h_ges+=(fraction*UnivGasConstant*t/comp.mol_wgt)*(-sp.cea_const[9]/t**2+sp.cea_const[10]*log(t)/t+sp.cea_const[11]+sp.cea_const[12]*t/2+sp.cea_const[13]*t**2/3+sp.cea_const[14]*t**3/4+sp.cea_const[15]*t**4/5+sp.cea_const[16]/t)
            return h_ges
    def tp2h(self,t,p,comp):
        return self.ideal(t,p,comp,'h')
        #return self.ideal(t,p,comp,'h')+self.dep(t,p,comp,'h)