'''
Created on 06.02.2017

@author: roysonntag
'''

import classes as cl


class Station():

    def __init__(self,t=None,p=None,w=None,comp=None):
        if t!=None and p!=None and w!=None and comp!=None:
            self.set_sation(t, p, w, comp)
        else:
            self.t=0
            self.p=0
            self.cp=0
            self.h=0
            self.s=0
            self.fluid = cl.CeaFluid()
            self.comp = cl.Composition()
            self.w=0 
    
    def calc_station(self):
        self.cp = self.fluid.tp2cp(self.t, self.p, self.comp)
        self.h = self.fluid.tp2h(self.t, self.p, self.comp)
        self.s = self.fluid.tp2s(self.t, self.p, self.comp)
        
    def set_sation(self,t,p,w,comp):
        self.t = t
        self.p = p
        self.w = w
        self.comp = comp