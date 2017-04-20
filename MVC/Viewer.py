import tkinter as tk
import tkinter.ttk as ttk
import matplotlib as mpl
mpl.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

from ctypes.test.test_pickling import name

class View(tk.Toplevel):
    def __init__(self, master=None):
        tk.Toplevel.__init__(self,master)
        # python dict holding widget handles
        self._widgets={}
        self._tree_col_id={}
        
        #### Window
        # destroy view with click on close button top left 
        self.protocol('WM_DELETE_WINDOW', self.master.destroy)
        self.title('RealGas')
        self.geometry('900x500')
        self.configure(background='gray')

        #### Axis
        self.figure = mpl.figure.Figure(figsize=(5,4), dpi=100)
        a = self.figure.add_subplot(111)
        a.plot([1,2,3,4,5,6,7,8],[5,6,1,3,8,9,3,5])
        self.canvas = FigureCanvasTkAgg(self.figure,self)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row=0,column=1,sticky='nsew',padx=1,pady=1)

        ### Initiate Note book
        self.add_widget('notebook', ttk.Notebook(self, height=400, width=400))
        self._widgets['notebook'].grid(row=0,column=0,sticky='nsew',padx=1,pady=1)
        
        ################## add tab for Fluid definition
        self.add_widget('fluid_tab',FluidTab(self._widgets['notebook'],'Fluids'))
        self._widgets['notebook'].add(self._widgets['fluid_tab'], text='Fluids', compound=tk.TOP)
        ################## add tab for Compression
        self.add_widget('compression_tab',CompressionTab(self._widgets['notebook'],'Compression'))
        self._widgets['notebook'].add(self._widgets['compression_tab'], text='Compression', compound=tk.TOP)
        ################## add tab for spec. heat
        self.add_widget('cp_tab',CpTab(self._widgets['notebook'],'CP'))
        self._widgets['notebook'].add(self._widgets['cp_tab'], text='CP', compound=tk.TOP)        
        
        ################## add fluid overview
        tree =  ttk.Treeview(self)    
        self.add_widget('table_composition',tree)
        # edit tree header
        self._widgets['table_composition'].heading('#0', text='Composition')
        # add standard species
        self.add_col_to_tree(['O2','CO2','H2O','AR','N2'])
        self._widgets['table_composition'].grid(row=1,columnspan=2,sticky='nsew',padx=1,pady=1)

    def add_widget(self,name, widget):
        self._widgets[name]=widget
        
    def add_composition(self,comp_list):
        # new comp is last entry of comp list
        new_comp = comp_list[-1]
        # add composition as new item
        id2 = self._widgets['table_composition'].insert("", 'end', "dir2", text=new_comp.name)
        # iterate through columns and set fraction for composition
        for species,col_id in self._tree_col_id.items():
            col = col_id
            fraction = new_comp.get_fraction(species)
            self._widgets['table_composition'].set(id2, column=col, value=fraction)
        #self._widgets['table_composition'].insert(id2, "end", "dir 2", text="test", values=("2A","2B"))
    
    def add_col_to_tree(self,names):
        col_tuple = self._widgets['table_composition']['columns']
        # check if tree is empty
        if len(col_tuple)==0:
            self._widgets['table_composition']['columns']=tuple(range(len(names)))
        else:
            self._widgets['table_composition']['columns']=col_tuple[:-1]+(col_tuple[-1],len(col_tuple))
        
        # iterate through names and start id at last id of tree
        for col_id, name in enumerate(names,len(col_tuple)):
            # set column width
            self._widgets['table_composition'].column(col_id,width=80)
            # set column header
            self._widgets['table_composition'].heading(col_id,text=name)
            # add new column to column dict
            self._tree_col_id[name]=col_id
    def update_figure(self,fig):
            self.figure = fig
            self.canvas = FigureCanvasTkAgg(self.figure,self)
            self.canvas.show()
        
            
class SimpleTable(tk.Frame):
    def __init__(self,parent,rows=5,columns=2):
        tk.Frame.__init__(self,parent,background='black')
        self._columns = columns
        self._rows = rows
        self._widgets=[]
        ### init headers
        current_row=[]
        for column in range(self._columns):
            label = tk.Label(self,text='head %d'%column,width=10,borderwidth=0)
            label.grid(row=0,column=column,sticky='nsew',padx=1,pady=1)
            current_row.append(label)
        self._widgets.append(current_row)

        ### init cells
        for row in range(1,self._rows):
            current_row=[]
            for column in range(self._columns):
                entry = tk.Entry(self,width=10,borderwidth=0)
                entry.grid(row=row,column=column,sticky='nsew',padx=1,pady=1)
                current_row.append(entry)
            self._widgets.append(current_row)

        for column in range(self._columns):
            self.grid_columnconfigure(column,weight=1)

    def set_cell(self,row,column,value):
        widget=self._widgets[row+1][column]
        widget.insert(0,value)
    def set_header(self,column,value):
        widget=self._widgets[0][column]
        widget.configure(text=value)
        
    def add_row(self, num_rows=1):
        self._rows+=num_rows
        for row in range(1,self._rows):
            current_row=[]
            for column in range(self._columns):
                entry = tk.Entry(self,width=10,borderwidth=0)
                entry.grid(row=row,column=column,sticky='nsew',padx=1,pady=1)
                current_row.append(entry)
            self._widgets.append(current_row)
        
        for column in range(self._columns):
            self.grid_columnconfigure(column,weight=1)
            
            

class Tab(tk.Frame):
    def __init__(self, master, name):
        tk.Frame.__init__(self, master)
        self.name = name
        self._rows=0
        self._widgets={}
        self.pack()
        self.bind("<Visibility>", self.on_visibility)

    def on_visibility(self, event):
        self.update()

    def add_widget(self,name, widget):
        self._widgets[name]=widget
        self._widgets[name].grid(row=self._rows, sticky='W')
        self._rows+=1
    
    def add_checkbox(self,name,text):
        self.add_widget(name, tk.Checkbutton(self, text=text))
        
    def add_entry(self,name,widget):
        self.add_widget(name,widget)
        
class FluidTab(Tab):
    def __init__(self, master, name):
        Tab.__init__(self, master, name)
        self._radio_var = {}
        ### radio buttons for fluid method
        self._radio_var['fluid'] = tk.IntVar()      
        self.add_widget('radiobutton_fluid_cea', tk.Radiobutton(self, text="CEA", variable=self._radio_var['fluid'], value=1))
        self.add_widget('radiobutton_fluid_rrd', tk.Radiobutton(self, text="RRD", variable=self._radio_var['fluid'], value=2))
        self.add_widget('radiobutton_fluid_bucker', tk.Radiobutton(self, text="Bucker", variable=self._radio_var['fluid'], value=3))      
        
        ### radio buttons for composition definition
        self._radio_var['comp'] = tk.IntVar()      
        self.add_widget('radiobutton_air', tk.Radiobutton(self, text="Air", variable=self._radio_var['comp'], value=1))
        self.add_widget('radiobutton_jet_a1', tk.Radiobutton(self, text="Jet-A1", variable=self._radio_var['comp'], value=2))
        self.add_widget('radiobutton_humid_air', tk.Radiobutton(self, text="Humid Air", variable=self._radio_var['comp'], value=3))
        self.add_widget('radiobutton_custom', tk.Radiobutton(self, text="Custom", variable=self._radio_var['comp'], value=4))
        
        self.add_entry('comp_name_entry',InputValue(self,text='Composition:'))

        comp_table = SimpleTable(self)
        comp_table.set_header(0,'Species')
        comp_table.set_header(1,'Vol. Fraction')
        self.add_widget('add_comp_table',comp_table)

        self.add_widget('add_comp_button',tk.Button(self, text='Add', width=8))
        self._widgets['add_comp_button'].grid(sticky='E')
        

        
class CompressionTab(Tab):
    def __init__(self, master, name):
        Tab.__init__(self, master, name)
        self.add_entry('t_in_entry',InputValue(self,text='T_in',unit='[K]'))
        self.add_entry('p_in_entry',InputValue(self,text='p_in',unit='[psi]'))
        self.add_entry('p_out_entry',InputValue(self,text='p_out',unit='[psi]'))

class CpTab(Tab):
    def __init__(self, master, name):
        Tab.__init__(self, master, name)

        self.add_entry('t_start',InputValue(self,text='T_start',unit='[K]',default=300))
        self.add_entry('t_end',InputValue(self,text='T_end',unit='[K]',default=2000))
        self.add_entry('p_in',InputValue(self,text='p_in',unit='[psi]',default=1e5))
        

class InputValue(tk.Frame):
    def __init__(self,parent,text,unit=None,default=None):
        tk.Frame.__init__(self,parent)
        if unit==None:
            self.label_name = tk.Label(parent,text=text)
            self.label_name.grid(row=parent._rows,column=0, sticky='W')
            
            self.entry = tk.Entry(parent,bd=5)
            self.entry.grid(row=parent._rows,column=1, sticky='W')
        else:
            self.label_name = tk.Label(parent,text=text)
            self.label_name.grid(row=parent._rows,column=0, sticky='W')
            
            self.entry = tk.Entry(parent,bd=5)
            self.entry.grid(row=parent._rows,column=1, sticky='W')
            
            self.label_unit = tk.Label(parent,text=text)
            self.label_unit.grid(row=parent._rows,column=2, sticky='W')
        if default:
            self.entry.insert(0,str(default))
        
  
        
        
        
    

