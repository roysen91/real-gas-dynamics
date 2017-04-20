import matplotlib as mpl
import matplotlib.pyplot as plt
#mpl.use('TkAgg')
import cantera as ct
import numpy as np
import classes as cl
import combust as cmb
import station as stat
from constants import FARStoi, AirToFuelFlow,PressureConversion, FLowConversion,StandartPressure,StandartTemperature
from sandbox import kappa_real



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

pressure_ratio = np.linspace(1,100,400)
#plot_work_press_ratio(pressure_ratio,StandartTemperature)
plot_eff_press_ratio(pressure_ratio,StandartTemperature)