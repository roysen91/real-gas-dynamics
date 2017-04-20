'''
Created on 03.02.2017

@author: roysonntag
'''

import unittest
import Controller as ctr
import classes as cl
from constants import *


class TestLoadCeaConstants(unittest.TestCase):
    def test_read_O2_cea_constants(self):
        my_species = cl.Species('O2')
        cea_const = [ -3.425563420E+04,4.847000970E+02,1.119010961E+00,4.293889240E-03,-6.836300520E-07,-2.023372700E-09,1.039040018E-12,-3.391454870E+03,1.849699470E+01,-1.037939022E+06,2.344830282E+03,1.819732036E+00,1.267847582E-03,-2.188067988E-07,2.053719572E-11,-8.193467050E-16,-1.689010929E+04,1.738716506E+01]
        self.assertEqual(cea_const,my_species.cea_const)
    def test_read_N2_cea_constants(self):
        my_species = cl.Species('N2')
        cea_const = [2.210371497E+04,-3.818461820E+02,6.082738360E+00,-8.530914410E-03,1.384646189E-05,-9.625793620E-09,2.519705809E-12,7.108460860E+02,-1.076003744E+01,5.877124060E+05,-2.239249073E+03,6.066949220E+00,-6.139685500E-04,1.491806679E-07,-1.923105485E-11,1.061954386E-15,1.283210415E+04,-1.586640027E+01]
        self.assertEqual(cea_const,my_species.cea_const)
    def test_read_CO2_cea_constants(self):
        my_species = cl.Species('CO2')
        cea_const = [4.943650540E+04,-6.264116010E+02,5.301725240E+00,2.503813816E-03,-2.127308728E-07,-7.689988780E-10,2.849677801E-13,-4.528198460E+04,-7.048279440E+00,1.176962419E+05,-1.788791477E+03,8.291523190E+00,-9.223156780E-05,4.863676880E-09,-1.891053312E-12,6.330036590E-16,-3.908350590E+04,-2.652669281E+01]
        self.assertEqual(cea_const,my_species.cea_const)
    def test_read_H2O_cea_constants(self):
        my_species = cl.Species('H2O')
        cea_const = [-3.947960830E+04,5.755731020E+02,9.317826530E-01,7.222712860E-03,-7.342557370E-06,4.955043490E-09,-1.336933246E-12,-3.303974310E+04,1.724205775E+01,1.034972096E+06,-2.412698562E+03,4.646110780E+00,2.291998307E-03,-6.836830480E-07,9.426468930E-11,-4.822380530E-15,-1.384286509E+04,-7.978148510E+00]
        self.assertEqual(cea_const,my_species.cea_const)
    def test_read_jeta1_L_cea_constants(self):
        my_species = cl.Species('Jet-A1(L)')
        cea_const = [ -4.218262130E+05,-5.576600450E+03, 1.522120958E+02, -8.610197550E-01,3.071662234E-03,-4.702789540E-06,2.743019833E-09,-3.238369150E+04,-6.781094910E+02]
        self.assertEqual(cea_const,my_species.cea_const)
    def test_read_jeta1_g_cea_constants(self):
        my_species = cl.Species('Jet-A1(g)')
        cea_const = [-6.068695590E+05 ,8.328259590E+03,-4.312321270E+01, 2.572390455E-01  ,  -2.629316040E-04 ,   1.644988940E-07  ,  -4.645335140E-11   , -7.606962760E+04  ,  2.794305937E+02  ,  1.858356102E+07 ,   -7.677219890E+04   , 1.419826133E+02  ,  -7.437524530E-03  ,  5.856202550E-07 ,   1.223955647E-11    ,-3.149201922E-15,    4.221989520E+05  ,  -8.986061040E+02]
        self.assertEqual(cea_const,my_species.cea_const)
        
class TestCreateComposition(unittest.TestCase):   
    def test_create_dry_air(self): 
        '''test if dry air is set up correctly'''
        my_comp = cl.Composition(comp_def='Air')
        self.assertEqual(my_comp.name, 'Air')
        self.assertEqual(len(my_comp.structure), 4)
        self.assertEqual(round(my_comp.r_spec,SignificantDigits), round(287.0689197920102,SignificantDigits))
        self.assertEqual(round(my_comp.mol_wgt,SignificantDigits), round(0.02896273134,SignificantDigits))
        self.assertEqual(my_comp.FAR, 0)
        self.assertEqual(my_comp.WGR, 0)
    def test_create_wet_air_with_RH(self): 
        '''test if wet air is set up correctly using relative humidity'''
        my_fluid = cl.CeaFluid()
        my_comp = cl.Composition(comp_def='Air')
        my_fluid.set_humidity(StandartTemperature+20, StandartPressure, my_comp, RH=1)
        self.assertEqual(my_comp.name, 'Air')
        self.assertEqual(len(my_comp.structure), 5)
        self.assertEqual(round(my_comp.r_spec,SignificantDigits), round(293.269688447603,SignificantDigits))
        self.assertEqual(round(my_comp.mol_wgt,SignificantDigits), round(0.028350355756201765,SignificantDigits))
        self.assertEqual(my_comp.FAR, 0)
        self.assertEqual(round(my_comp.WGR,SignificantDigits), round(0.036854216859696774,SignificantDigits))
        self.assertEqual(round(my_comp.get_fraction('H2O'),SignificantDigits),  round(0.05593629993099881,SignificantDigits))
    def test_create_wet_air_with_WGR(self): 
        '''test if wet air is set up correctly using WGR'''
        my_fluid = cl.CeaFluid()
        my_comp = cl.Composition(comp_def='Air')
        my_fluid.set_humidity(StandartTemperature+20, StandartPressure, my_comp, WGR=0.036854216859696774)
        self.assertEqual(my_comp.name, 'Air')
        self.assertEqual(len(my_comp.structure), 5)
        self.assertEqual(round(my_comp.r_spec,SignificantDigits), round(293.269688447652,SignificantDigits))
        self.assertEqual(round(my_comp.mol_wgt,SignificantDigits), round(0.028350355756201765,SignificantDigits))
        self.assertEqual(my_comp.FAR, 0)
        self.assertEqual(round(my_comp.WGR,SignificantDigits), round(0.036854216859696774,SignificantDigits))
        self.assertEqual(round(my_comp.get_fraction('H2O'),SignificantDigits),  round(0.05593629993099881,SignificantDigits))
    def test_create_jeta1(self): 
        '''test if Jet-A1 is set up correctly'''
        my_comp = cl.Composition(comp_def='Jet-A1')
        self.assertEqual(my_comp.name, 'Jet-A1')
        self.assertEqual(len(my_comp.structure), 1)
        self.assertEqual(round(my_comp.r_spec,SignificantDigits), round(48.878894767783656,SignificantDigits))
        self.assertEqual(my_comp.mol_wgt, 0.1701)
        self.assertEqual(my_comp.FAR, 1)
        self.assertEqual(my_comp.WGR, 0)
        

            
#class TestAddComposition(unittest.TestCase):
    
    #def test_air(self):
    #    root = ctr.tk.Tk()
    #    root.withdraw()
    #    app=ctr.Controller(root)
    #    app.view._widgets['fluid_tab']._radio_var.set(1)
    #    app.add_composition()
    #   self.assertEqual(app.view._widgets['table_composition'].item('end'), 'Air')

#class TestIndependent(unittest.TestCase):

 #   def test_empty(self):
  #      X = []
   #     e = [1, 0]
    #    self.assertTrue(is_independent_of(e, X))
     #   self.assertTrue(e in X)

if __name__ == '__main__':
    app=ctr.Controller()
    unittest.main()
