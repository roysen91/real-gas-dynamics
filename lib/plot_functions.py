# Import all necessary libraries
import matplotlib as mpl
import matplotlib.pyplot as plt
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
        if comp.WGR_mole>0:
            file1='data/CEA_tool/cea_full_rh=1.dat'
            with open(file1,'r') as file1_handle:
                # skip first row
                next(file1_handle)
                for line in file1_handle:
                    line=line.strip().split('  ')
                    eqr_range_cea.append(float(line[0]))
                    t_onl_cea.append(float(line[1]))
                file1_handle.close()   
        elif comp.WGR_mole==0:
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
        diff_t_rrd = (t_rrd-t_base)/t_base*100
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
        diff_t_cea = (t_cea-t_base)/t_base*100
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
        diff_t_cantera = (t_cantera-t_base)/t_base*100
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

def plot_cp(t_range = np.linspace(300,2000,100),p=1e5,RH = 0):
    fl_cea = cl.CeaFluid()
    fl_rrd = cl.RRDFluid()
    std_air= cl.Composition('Air')

    fl_cea.set_humidity(StandartTemperature+20, StandartPressure, std_air, RH = RH)

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

def plot_eff_press_ratio(pressure_ratio,t_in):
    comp = cl.Composition()
    fluid_cea = cl.CeaFluid(eos_model='SOA')
    fluid_cea.eos.set_coeffs(comp)
    eta_is = 0.9
    
    eta_const = []
    eta_real = []
    t_out_test =[]
    t_out_ref = []
    t_out_real = []
    eta_base1=[]
    eta_base2=[]
    
    for pr in pressure_ratio:
        #if pr == 1:
        #    eta_const.append(eta_is)
        #    eta_real.append(eta_is)
        #else:
        h_in = fluid_cea.tp2h(t_in,StandartPressure,comp)
        kappa = fluid_cea.tp2kappa(t_in, StandartPressure, comp)
        t_out_test.append(t_in*(((pr**((kappa-1)/kappa)-1)/eta_is)+1))
        t_out_ref.append(fluid_cea.isentropic_compression(t_in, StandartPressure, pr*StandartPressure, comp, eta_isentropic=eta_is))
        t_out_real.append(fluid_cea.isentropic_compression(t_in, StandartPressure, pr*StandartPressure, comp,mode='real', eta_isentropic=eta_is))
        
        kappa_mean = (kappa+fluid_cea.tp2kappa(t_out_ref[-1],pr*StandartPressure,comp))/2
        
        eta_base1.append((pr**((kappa-1)/kappa)-1)/((t_out_test[-1]/t_in)-1))
        eta_base2.append((pr**((kappa_mean-1)/kappa_mean)-1)/((t_out_ref[-1]/t_in)-1))
        
        s_is = fluid_cea.tp2s(t_in,StandartPressure,comp)
        t_out_ref_is = fluid_cea.sp2t(s_is,pr*StandartPressure,comp)
        h_out_ref_is = fluid_cea.tp2h(t_out_ref_is,pr*StandartPressure,comp)
        #h_out_const = fluid_cea.tp2h(t_out_ref[-1],pr*StandartPressure,comp)
        #h_out_real = fluid_cea.tp2h(t_out_ref[-1],pr*StandartPressure,comp,mode='real')
        h_out_const = fluid_cea.tp2h(t_out_real[-1],pr*StandartPressure,comp)
        h_out_real = fluid_cea.tp2h(t_out_real[-1],pr*StandartPressure,comp,mode='real')
        

        eta_const.append((h_out_ref_is-h_in)/(h_out_const-h_in))
        eta_real.append((h_out_ref_is-h_in)/(h_out_real-h_in))
        #eta_const.append((t_out_ref_is-t_in)/(t_out_real[-1]-t_in))
    eta_const=np.array(eta_const)  
    eta_real=np.array(eta_real)  
    t_out_test=np.array(t_out_test) 
    t_out_ref=np.array(t_out_ref) 
    eta_base1=np.array(eta_base1)
    eta_base2=np.array(eta_base2)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    #ax1.plot(pressure_ratio,eta_const,pressure_ratio,eta_mean_ideal,pressure_ratio,eta_mean_real)
    ax1.plot(pressure_ratio,eta_const,pressure_ratio,eta_real,[pressure_ratio[0],pressure_ratio[len(pressure_ratio)-1]],[0.9,0.9])
    ax1.set_xlabel('Pressure ratio [1]')
    ax1.set_ylabel('isentropic efficiency [1]')
    ax1.grid(True)
    ax1.legend(['ideal $\eta_{is}$','real $\eta_{is}$'])


    '''
    ax2 = fig.add_subplot(2,2,2)
    ax2.plot(pressure_ratio,(eta_real-eta_const))
    ax2.set_xlabel('Pressure ratio [1]')
    ax2.set_ylabel('delta isentropic efficiency [1]')
    ax2.grid(True)
    ax2.legend(['dev. real $\eta_{is}$'])
    
    ax3 = fig.add_subplot(2,2,3)
    ax3.plot(pressure_ratio,t_out_test/t_in,pressure_ratio,t_out_ref/t_in)
    ax3.set_xlabel('Pressure ratio [1]')
    ax3.set_ylabel('Outlet Temperature [K]')
    ax3.grid(True)
    ax3.legend(['$\kappa=const$','$\kappa=(\kappa_{in}+\kappa_{out})/2$'])
    
    ax4 = fig.add_subplot(2,2,4)
    ax4.plot(pressure_ratio,eta_base1,pressure_ratio,eta_base2)
    ax4.set_xlabel('Pressure ratio [1]')
    ax4.set_ylabel('isentropic efficiency [1]')
    ax4.set_ylim([0.85,0.95])
    ax4.grid(True)
    ax4.legend(['$\kappa=const$','$\kappa=(\kappa_{in}+\kappa_{out})/2$'])
    '''
    plt.show()
    
    #plot_t_out(t_in,t_out_test,t_out_ref,t_out_real)
    plot_pr_kappa(t_in,pressure_ratio)
    
def plot_pr_kappa(t_in,pressure_ratio):
    comp = cl.Composition()
    fluid_cea = cl.CeaFluid(eos_model='SOA')
    fluid_cea.eos.set_coeffs(comp)
    eta_is = 0.9
    
    kappa_mean =[]
    kappa_ideal = fluid_cea.tp2kappa(t_in, StandartPressure, comp)
    
    f1 = []
    f2 = []
    
    for pr in pressure_ratio:
        t_out = fluid_cea.isentropic_compression(t_in, StandartPressure, pr*StandartPressure, comp, eta_isentropic=eta_is)
        kappa_mean.append((kappa_ideal+fluid_cea.tp2kappa(t_out,pr*StandartPressure,comp))/2)
        f1.append(pr**((kappa_ideal-1)/kappa_ideal))
        f2.append(pr**((kappa_mean[-1]-1)/kappa_mean[-1]))
        
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    
    ax1.plot(pressure_ratio,f1,pressure_ratio,f2)
    ax1.set_xlabel('Pressure ratio [1]')
    ax1.set_ylabel('f(pressure ratio)')
    #ax1.set_xticks([1,20,30,40,50,60,70,80,90,100])
    #ax1.set_yticks([StandartTemperature,400,600,800,1000,1200,1400])
    ax1.grid(True)
    ax1.legend(['$\kappa_{ideal}=const$','$\kappa_{ideal}=(\kappa_{in}+\kappa_{out})/2$'])
    
    plt.show()
        
def plot_t_out(t_in,t_out_test,t_out_ref,t_out_real):
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    
    ax1.plot(pressure_ratio,t_out_test,pressure_ratio,t_out_ref,pressure_ratio,t_out_real)
    ax1.set_xlabel('Pressure ratio [1]')
    ax1.set_ylabel('Outlet Temperature [K]')
    #ax1.set_xticks([1,20,30,40,50,60,70,80,90,100])
    #ax1.set_yticks([StandartTemperature,400,600,800,1000,1200,1400])
    ax1.grid(True)
    ax1.legend(['$\kappa_{ideal}=const$','$\kappa_{ideal}=(\kappa_{in}+\kappa_{out})/2$','$\kappa_{real}=(\kappa_{in}+\kappa_{out})/2$'])
    
    plt.show()
        
def plot_work_press_ratio(pressure_ratio,t_in):
    comp = cl.Composition()
    fluid_cea = cl.CeaFluid(eos_model='SOA')
    fluid_cea.eos.set_coeffs(comp)
    fluid_rrd = cl.RRDFluid()
    
    delta_h_const = []
    delta_h_mean_ideal = []
    delta_h_mean_real = []
    
    for pr in pressure_ratio:
        kappa = fluid_cea.tp2kappa(t_in, StandartPressure, comp)
        t_out_const = t_in *(pr)**((kappa-1)/kappa)
        t_out_mean_ideal = fluid_cea.isentropic_compression(t_in, StandartPressure, pr*StandartPressure, comp)
        t_out_mean_real = fluid_cea.isentropic_compression(t_in, StandartPressure, pr*StandartPressure, comp,mode='real')
        delta_h_const.append(fluid_cea.tp2h(t_out_const, StandartPressure, comp)-fluid_cea.tp2h(t_in, pr*StandartPressure, comp))
        delta_h_mean_ideal.append(fluid_cea.tp2h(t_out_mean_ideal, StandartPressure, comp)-fluid_cea.tp2h(t_in, pr*StandartPressure, comp))
        delta_h_mean_real.append(fluid_cea.tp2h(t_out_mean_real, StandartPressure, comp,mode='real')-fluid_cea.tp2h(t_in, pr*StandartPressure, comp,mode='real'))
    
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    #ax1.plot(pressure_ratio,eta_const,pressure_ratio,eta_mean_ideal,pressure_ratio,eta_mean_real)
    ax1.plot(pressure_ratio,delta_h_const,pressure_ratio,delta_h_mean_ideal,pressure_ratio,delta_h_mean_real)
    ax1.set_xlabel('Pressure ratio [1]')
    ax1.set_ylabel('spec. compressor work [J/kg]')
    ax1.grid(True)
    ax1.legend(['$\delta h$ mit $\kappa=const$','$\delta h_{id}$ mit $\kappa=(T_{in}+T_{out})/2$','$\delta h_{re}$ mit $\kappa=(T_{in}+T_{out})/2$'])
    plt.show()
def plot_burn_temp_comparison(p_in,t_in,w_in,t_fuel=382.5948):
    # load large and reduced chemestry set data from CEA website
    FileNames=['data/online_cea_combustion/reducedChemSet_EQR_below1.dat',
               'data/online_cea_combustion/reducedChemSet_EQR_above1.dat',
               'data/online_cea_combustion/largeChemSet_EQR_below1.dat',
               'data/online_cea_combustion/largeChemSet_EQR_above1.dat']

    # read data from loaded files
    T = np.array([])
    CP = np.array([])
    EQR = np.array([])
        
    for i,inFile in enumerate(FileNames):
        if i==2:
            Tred=T
            CPred=CP
            EQRred=EQR
            
            T = np.array([])
            CP = np.array([])
            EQR = np.array([])
        inputFile = open(inFile,'r')

        for line in inputFile:
            if 'T, K' in line:
                line = line.split()
                T = np.append(T,float(line[2]))
            if 'Cp, KJ/(KG)(K)' in line:
                line = line.split()
                CP = np.append(CP,float(line[2]))
            if 'EQ.RATIO' in line:
                line = line.split(',')
                EQR = np.append(EQR,float(line[-1].split('=')[1]))

        inputFile.close()

    # reduced set of species using my cea calculation
    dummy_fluid = cl.CeaFluid()
    std_air=cl.Composition('Air')
    dummy_fluid.set_humidity(StandartTemperature+20, StandartPressure, std_air, RH=1)

    # set up combustor
    comb = cmb.Combustor()
    comb.current_fluid = cl.CeaFluid()
    comb.in_station = stat.Station(t_in,p_in*PressureConversion,w_in*FLowConversion,std_air)
    comb.out_station = stat.Station(t_in,p_in*PressureConversion,w_in*FLowConversion,std_air)
    comb.current_fuel.tf = t_fuel

    eqr_vec = np.linspace(0.25, 2, 20)
    far_vec = eqr_vec*FARStoi
    wf_vec = far_vec*comb.in_station.w*AirToFuelFlow

    t_cea = np.array([])

    for fuel_flow in np.nditer(wf_vec):
        comb.wf = fuel_flow
        comb.heat_balance_calcs()
        t_cea = np.append(t_cea,comb.out_station.t)
        
    comb.current_fluid = cl.RRDFluid()
    comb.current_fuel.burn_model = 'RRD'

    t_rrd = np.array([])

    for fuel_flow in np.nditer(wf_vec):
        comb.wf = fuel_flow
        comb.heat_balance_calcs()
        t_rrd = np.append(t_rrd,comb.out_station.t)

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
        gas2.TP = t_in, p_in*PressureConversion
        gas2.set_equivalence_ratio(EQR[i], 'Jet-A(g)', 'O2:0.19775868763, N2:0.7371626995, AR:0.008841156551, CO2:0.0003011563203 , H2O: 0.05593629993')
        gas2.equilibrate('HP')
        T_incomplete[i] = gas2.T

    ##### PLOT 
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