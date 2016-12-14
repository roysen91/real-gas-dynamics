'''
Created on 17.11.2016

@author: roysonntag
'''
from classes import *
from Viewer import *

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

class Controller:
	def __init__(self):
		root = tk.Tk()
		root.withdraw()

		self.view = View(root)
		self.view.mainloop()

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

Controller()