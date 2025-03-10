#!/usr/bin/env python3 

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import csv

'''
Lattice_Plot.py version 0.0.2

Creates Lattice_Plots as spacial maps of order parameters. 
To be used with data created by "Lipid_Lattice.c".
Namely output file "lattice.csv".

Example Input and Output Files can be found in "/Check" for comparison after
changing the code.

Changelog:

0.0.2:  Changed input from ".dat" to ".csv"
        Added example files in "/Check"

0.0.1:  Initial release.
'''

class Lattice_Plot:
    '''Takes lattice data and creates a spacial map of order parameters and lipid types.'''
    def __init__(self,dim:int = 0.0) -> None:
        '''Initializes a lattice with dimension dim.'''
        self.dimension = dim
        self.lattice_order = np.zeros(shape=(dim,dim))
        self.lattice_types = np.zeros(shape=(dim,dim))
        self.lattice_chol = []
        pass

    def read_file(self,path:str) -> None:
        '''Reads the data from a lattice.csv file.'''  

        with open(path) as file:
            lines = [line.rstrip() for line in file]

        vector = [0 for x in range(self.dimension**2)]

        for i in range(self.dimension**2):
            vector[i] = lines[i].split(",")

        for coord in vector:
            x = int(coord[0])
            y = int(coord[1])
            t = int(coord[2])
            s = int(coord[3])
            self.lattice_types[x][y] = t
            self.lattice_order[x][y] = (s-50)/100

        vector_chol = [0 for x in range(len(lines) - self.dimension**2)]
        x_chol = []
        y_chol = []

        for i in range(len(lines)-self.dimension**2):
            vector_chol[i] = lines[i + self.dimension**2].split(",")

        for coord in vector_chol:
            x_chol.append(float(coord[0]))
            y_chol.append(float(coord[1]))

        file.close()

        self.lattice_chol = [x_chol,y_chol]

    def plot_types(self,path:str="Lipid_Type.png") -> None:
        '''Plots the spacial lattice with respect to different lattice types.'''

        x = self.lattice_chol[1]
        y = self.lattice_chol[0]

        fig,ax = plt.subplots()

        cmap = colors.ListedColormap(["cornflowerblue",
                                      "lightgreen"])
        bounds = [1,2,3]
        norm = colors.BoundaryNorm(bounds,cmap.N)

        im = ax.imshow(self.lattice_types,cmap=cmap,norm=norm)
        ax.scatter(x,y,marker="o",s=1,c="black",edgecolors="black")
        ax.grid(False)

        plt.colorbar(im,ax=ax)
        plt.savefig(path)

    def plot_orders(self,path:str="Lipid_Order.png") -> None:
        '''Plots the spacial lattice with respect to the lipids' order parameter.'''

        x = self.lattice_chol[1]
        y = self.lattice_chol[0]

        color_lattice = self.map_colors(self.lattice_order)

        fig,ax = plt.subplots()

        im = ax.imshow(color_lattice,vmin=-0.5,vmax=1.0)
        ax.scatter(x,y,marker="o",s=1,c="black",edgecolors="black")
        ax.grid(False)

        cbar = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.05)
        cbar.ax.invert_xaxis()
        cbar.set_ticks([-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])
        cbar.set_ticklabels([-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1])

        plt.savefig(path)


    def map_colors(self,arr:list) -> list:
        '''Generates a colormap for the lipids' order parameter.'''
        colored_array = np.empty((*arr.shape, 3))
        norm = colors.Normalize(vmin=-0.5,vmax=1.0)
        colormap = plt.get_cmap("viridis")

        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                value = arr[i,j]
                colored_array[i,j] = colormap(norm(value))[:-1]

        return colored_array

def main(files:list,names:list):
    '''Plots the lattices of a given file and label list.'''
    lattice = Lattice_Plot(10)
    for index,file in enumerate(files):
        lattice.read_file(file)
        lattice.plot_orders(names[index]+"_Order.png")
        lattice.plot_types(names[index]+"_Types.png")

files = ["0/lattice.csv"]

names = ["Test"]

main(files,names)