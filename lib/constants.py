'''
Created on 06.02.2017

@author: roysonntag
'''
UnivGasConstant= 8.3143 # [J/mol*K]
StandartPressure = 101325 # [Pa]
StandartTemperature = 288.15 # [K]
RRDSpecHeatConversion = 1.643989e+03 #[HPs/lb*K] to [J/kg*K]
FARStoi = 0.06810617
AirToFuelFlow = 3600 # [1]
MinTemp = 150 # [K]
MaxTemp = 3000 # [K]
MinMolVol = 5e-5
SpecHeatConversion = 1 # 0.000238805 [Chu/lbm*K] to [J/kg*K]
PressureConversion = 6894.757 # [psi] to [Pa]
FLowConversion = 0.453592 # [lbm/s] to [kg/s]
SignificantDigits = 8 # for rounding float numbers 
MaxIterations = 40
MinDeviation = 1e-7
StartValueMolareVolume = 5e-3