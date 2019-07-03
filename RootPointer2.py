#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 14:44:06 2019

@author: erich

Prescription: .VARIABLES, .MethodsUsed()

"""

import numpy as np
class lattice:
    def __init__(self,LL):
        """ 
        Initialize an empty lattice.
        """
        self.L = LL
        # lattice size L*L
        self.GRID = np.zeros((LL,LL),dtype=bool)
        # grid of occupied cells
        self.DIRX = np.zeros((LL,LL),dtype=int)
        self.DIRY = np.zeros((LL,LL),dtype=int)
        # grids that give the pointers to the root
        self.ListIndex()
        # generates a list of cells
        self.ROOTSTATS = {}
        # dictionary of roots and statistics
        
        self.POSX,self.POSY = np.mgrid[0:LL,0:LL]
        #static grid useful for .DefRootOfCell
        return
    def ListIndex(self):
        """
        Get a randomized list of cells. This is used to
        increase the occupacy number at each step
        """
        self.ALLINDEX = list(zip(*np.where(self.GRID==0)))
        np.random.shuffle(self.ALLINDEX)
        return
    def AddCell(self):
        LL = self.L
        "Adds a new cell to the grid using .ALLINDEX"        
        if len(self.ALLINDEX)==1: return
        # check if .ALLINDEX is empty
        INDEX = self.ALLINDEX.pop(0)
        self.GRID[INDEX] = True 
        # get one index, active the respective cell
    
        "Check roots of active neighborhood"        
        LIST_NBH_ROOT = self.ListNeighborRoot(INDEX)
        # gets the list of roots for each cell in the neighborhood
        LIST_NBH_ROOT, REPEATED_ROOT = self.CheckDuplicated(LIST_NBH_ROOT)
        # checks for duplicates 
        "Define the pointer of the cell using .AttachTo(). If it points to itself DIR(X|Y) = 0"
        LEN = len(LIST_NBH_ROOT)
        if LEN==0 : self.ROOTSTATS[INDEX]=[1,False] #Isolated cluster
        if LEN >0:
            self.DefRootOfCell(LIST_NBH_ROOT[0],INDEX)
            #Points to the first option
            ROOT_1st = (LIST_NBH_ROOT[0][0]%LL,LIST_NBH_ROOT[0][1]%LL)
            self.ROOTSTATS[ROOT_1st][0]+=1
            #Adds cell to a cluster
            ### MUST modify to POINT to the BIGGER Cluster
            for NGH_ROOT in LIST_NBH_ROOT[1:]:
                self.DefRootOfCell(LIST_NBH_ROOT[0],NGH_ROOT)
                #It makes all other roots point to the new root
                aux = self.ROOTSTATS.pop((NGH_ROOT[0]%LL,NGH_ROOT[1]%LL)) 
                #remove from stats
                self.ROOTSTATS[ROOT_1st][0]+=aux[0]
                self.ROOTSTATS[ROOT_1st][1]+=aux[1]
                #adds info to the stats of the new cluster
        if LEN>1 : self.RefreshDIR()
        #refresh grid IF we have a union between two different clusters
        ### MIGHT be implemented so that we do not refresh the WHOLE grid

        """Check percolation
        self.ROOTSTATS guarda informação se percolou 
        {ROOT INDEX at first sheet : [Cluster size, Percolated cluster]}
        {int tuple : [int, bool] }
        
        When adding a new cell one has some possibilities
        - At least one of the neighbors have percolated
            Therefore we must simply unite all these clusters and identify 
            that they have percolated. This is already done in former.
        - None has percolated and all roots are different.
            Therefore no "percolation check" is needed and the union is
            already done in former.
        - None has percolated BUT we have repeated roots.
            There is a need to check for percolation
        """
        if len(REPEATED_ROOT)>0:  #at least one repeated root
            if (self.ROOTSTATS[ROOT_1st][1]==False): # it has not percolated yet
                #note that ROOT_1st IS defined if LEN>0, which is guarantee
                #by len(REPEATED_ROOT)>0 
                if self.CheckPercolation(INDEX): #True if percolated
                    self.ROOTSTATS[ROOT_1st][1] = True
        return
    def ListNeighborRoot(self,INDEX):
        """
        List the root of the neighborhood of the given cell
        
        First it defines the neighbor index with AUX. 
        Then, if AUX is occupied it checks its root.
        
        To diminish the amount of times BC_P is called it is defined BAUX.
        """
        NB_R = [] 
        GRID0 = self.GRID
    
        AUX = (INDEX[0]-1,INDEX[1])
        BAUX = self.BC_P(AUX)
        if GRID0[BAUX] : NB_R += [self.RootOf(AUX,BAUX)]

        AUX = (INDEX[0]+1,INDEX[1])             
        BAUX = self.BC_P(AUX)
        if GRID0[BAUX] : NB_R += [self.RootOf(AUX,BAUX)]

        AUX = (INDEX[0],INDEX[1]-1) 
        BAUX = self.BC_P(AUX)
        if GRID0[BAUX] : NB_R += [self.RootOf(AUX,BAUX)]

        AUX = (INDEX[0],INDEX[1]+1) 
        BAUX = self.BC_P(AUX)
        if GRID0[BAUX] : NB_R += [self.RootOf(AUX,BAUX)]

        return NB_R
    def RootOf(self,INDEX,BINDEX):
        """
        Says what is the root of a given cell based on the sum of its
        location and its pointers.
        """
        return (INDEX[0]+self.DIRX[BINDEX],INDEX[1]+self.DIRY[BINDEX])
    def BC_P(self,INDEX):
        """
        Imposed boundary condition. In this case a periodic one.
        """
        LL = self.L
        return (INDEX[0]%LL,INDEX[1]%LL)
    def CheckDuplicated(self,LIST):
        """
        Takes a list of roots and check which are the uniques and which are
        repeated.
        To do so it MUST apply boundary conditions of the index to see if they
        are equivalent
        
        Maybe the code can be incremented using set(). 
        But we have just 4 neighbors at max, so it seems unnecessary
        https://stackoverflow.com/questions/9835762/how-do-i-find-the-duplicates-in-a-list-and-create-another-list-with-them
        """
        LL = self.L
        UNIQUE = LIST
        REPEATED = []
        
        if len(LIST)>0:
            #only works if there is a list
            BC_LIST = [(ELEM[0]%LL,ELEM[1]%LL) for ELEM in LIST]
            #roots are now at the "main" sheet and can be compared
            for ELEM in BC_LIST:
                while BC_LIST.count(ELEM)>1 :
                    UNIQUE.pop( BC_LIST.index(ELEM) )
                    BC_LIST.remove(ELEM)
                    if ELEM not in REPEATED : REPEATED += [ELEM]
                    
        return UNIQUE,REPEATED #final information IS NOT using the main sheet
    def DefRootOfCell(self,ROOT,INDEX):
        "Creates a pointer to the root"
        BINDEX = self.BC_P(INDEX)
        self.DIRX[BINDEX] = ROOT[0]-INDEX[0]
        self.DIRY[BINDEX] = ROOT[1]-INDEX[1]
        return
    def RefreshDIR(self):
        "Refresh matrix os pointers."
        "Other option would be a 'recurrence' of DefRootOfCell to AVOID matrix operations"
        LL = self.L

        RootX = ((self.DIRX + self.POSX)*self.GRID)%LL
        RootY = ((self.DIRY + self.POSY)*self.GRID)%LL
        #Matrix of fake roots (first approximation to the root)
        BASE = self.GRID*((self.DIRX!=0)+(self.DIRY!=0)).any()
        self.DIRX += BASE*self.DIRX[(RootX,RootY)]
        self.DIRY += BASE*self.DIRY[(RootX,RootY)]
        #Pointer to the real root        
        return
    def CheckPercolation(self,INDEX):
        # check if difference of distance for any of the first neighboors is greater than 1.
        # This difference signals the path difference to the root and, therefore, the existence of percolation
        AUX = [(INDEX[0]-1,INDEX[1]),(INDEX[0]+1,INDEX[1]),(INDEX[0],INDEX[1]-1),(INDEX[0],INDEX[1]+1)]
        DIST = []
        DX0 = self.DIRX[INDEX]
        DY0 = self.DIRY[INDEX]
        for EACH in AUX:
            BEACH = self.BC_P(EACH)
            if self.GRID[BEACH]:
                DX = self.DIRX[BEACH] - DX0
                DY = self.DIRY[BEACH] - DY0
                DIST += [DX*DX+DY*DY]
        return any(np.array(DIST)>1.)
    
from matplotlib import pyplot as plt
def Vplot(grid,vecx,vecy):
    plt.figure()
    plt.imshow(grid,vmin=0,vmax=2,cmap='coolwarm')
    plt.quiver(vecy,-vecx,color='green',units='xy')
    plt.axis('off')
    return

###########################################################
##Simple test of code's well behavior
#Q = lattice(5)
#Q.AddCell()
#Vplot(Q.GRID*1 + Q.GRID*(Q.DIRX==0)*(Q.DIRY==0)*1 ,Q.DIRX,Q.DIRY)
#Q.ROOTSTATS

def one_run(L): 
    """"One sample"""
    N = L*L
    Q = lattice(L)
    STATS = []
    
    for n in range(N-1): #modifying the number of occupied cells
        Q.AddCell()
        STATS.append(Q.ROOTSTATS)

    return STATS

###########################################################
## Testing elapsed time
#import time
#
#L_range = np.linspace(5,300,20,dtype=int)
#T_range = np.zeros_like(L_range,dtype=float)
#
##for i in range(20):
#i=13
#T_range[i] = time.time()
#one_run(L_range[i])
#T_range[i] = time.time()-T_range[i]
#
#def powlaw(x,a,b):
#    return a*np.power(x,b)
#
#x = L_range[0:i+1]
#y = T_range[0:i+1]
#from scipy.optimize import curve_fit
#popt,pcov = curve_fit(powlaw,x,y)
#
#xlin = np.linspace(0,x[i],100)
#fig = plt.figure()
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
#ax.plot(x, y,'o:')
#ax.plot(xlin,powlaw(xlin,*(popt)),'r-',label='%5.1E $x^{%5.2f}$'% tuple(popt))
#ax.legend(loc="upper left")
#ax.set_xlabel('Grid size (L)')
#ax.set_ylabel('Run time (s)')
