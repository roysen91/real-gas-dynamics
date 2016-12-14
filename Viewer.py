import tkinter as tk

class View(tk.Toplevel):
    def __init__(self, master=None):
    	tk.Toplevel.__init__(self,master)
    	
    	#### Window
    	# destroy view with click on close button top left 
    	self.protocol('WM_DELETE_WINDOW', self.master.destroy)
    	self.title('RealGas')
    	self.geometry('1000x500')
    	self.configure(background='gray')

    	#### Axis
    	axis = tk.Canvas(self,bg="blue", height=250, width=300)
    	axis.pack(side=tk.RIGHT)
