'''
Created on 17.11.2016

@author: roysonntag
'''

from Viewer import *
from Model import *

class Controller:
	def __init__(self,root):
		self.model = Model()
		# add functions to observable
		self.model.fluids.add_callback(self.fluid_changed)
		self.model.fluids.add_callback(self.tab_changed)
		
		self.view = View(root)
		self.view._widgets['fluid_tab']._widgets['add_comp_button'].config(command=self.add_fluid)
		self.view._widgets['cp_tab'].config(command=self.tab_changed(self.view._widgets['cp_tab'].name))
		
		
	# gets called on click of Add button within Fluid Tab
	def add_fluid(self):
		fluid_tab_dict = {}
		fluid_tab_dict['type']=self.get_radio_button_value('fluid')
		fluid_tab_dict['comp']=self.get_radio_button_value('comp')
		self.model.add_fluid(fluid_tab_dict)

	# checks if a radio button is ticked	
	def radio_button_checked(self,rb_group):		
		if self.view._widgets['fluid_tab']._radio_var[rb_group].get() != 0:
			return True
		else:
			return False
	def get_radio_button_value(self,rb_group):
		if rb_group == 'fluid':
			if self.view._widgets['fluid_tab']._radio_var[rb_group].get() == 1:
				return 'CEA'
			elif self.view._widgets['fluid_tab']._radio_var[rb_group].get() == 2:
				return 'RRD'
			elif self.view._widgets['fluid_tab']._radio_var[rb_group].get() == 3:
				'Bucker'
		elif rb_group == 'comp':
			if self.view._widgets['fluid_tab']._radio_var[rb_group].get() == 1:
				return 'Air'
			elif self.view._widgets['fluid_tab']._radio_var[rb_group].get() == 2:
				return 'Jet-A1'

		
	# updates view if list of working fluids changes	
	def fluid_changed(self,comp):
		self.view.add_composition(comp)
	def tab_changed(self,tab_name):
		if tab_name == 'cp_tab':
			self.view.update_figure(self.model.plot_cp())
	
		
	
	#def plot_cp(x_axis_list):

		#f = Figure(figsize=(5,5), dpi=100)
		#a = f.add_subplot(111)
		#a.plot([1,2,3,4,5,6,7,8],[5,6,1,3,8,9,3,5])

		#canvas = FigureCanvasTkAgg(f, self)
		#canvas.show()
		#canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

		#toolbar = NavigationToolbar2TkAgg(canvas, self)
		#toolbar.update()
		#canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

# only gets executed if file is run by itself  and not executed
if __name__ == '__main__':
	root = tk.Tk()
	root.withdraw()
	app = Controller(root)
	root.mainloop()