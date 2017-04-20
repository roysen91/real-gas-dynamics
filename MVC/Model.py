import classes as cl
import matplotlib as mpl
import matplotlib.pyplot as plt

class Model():
	def __init__(self):
		# list holding fluids used for calculations
		self.fluids = Observable([])
	def add_fluid(self,fluid_dict):
		self.fluids.add(cl.Composition('Mixture_'+str(len(self.fluids.data)),comp_def))
	def plot_cp(self,t_range,p,compositions):
		fig = plt.figure()
		fluid = cl.CeaFluid()
		ax = fig.add_subplot(1,2,1)
		cp={}
		for comp in compositions:
			for temp in t_range:
				cp[comp]=fluid.tp2cp(temp,p,comp)
			ax.plot(t_range,cp[comp])
		return fig

		
class Observable():
	def __init__(self, initialValue=None):
		self.data = initialValue
		self.callbacks = {}
		
	def add_callback(self, func):
		self.callbacks[func] = 1

	def del_callback(self, func):
		del self.callback[func]

	def _do_callbacks(self):
		for func in self.callbacks:
			func(self.data)

	def set(self, data):
		self.data = data
		self._do_callbacks()
	
	def add(self,data):
		self.data.append(data)
		self._do_callbacks()

	def get(self):
		return self.data

	def unset(self):
		self.data = None
