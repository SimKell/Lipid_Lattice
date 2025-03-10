#!/usr/bin/env python3 

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
import seaborn as sns
import csv 

'''
Kite_Plot.py version 0.0.2

Creates Kite_Plots of the given order parameter distributions. 
To be used with data created by "Lipid_Lattice.c".
Namely output files "order.csv", "order_neigh.csv" and "order_corr.csv".

Example Input and Output Files can be found in "/Check" for comparison after
changing the code.

Changelog:

0.0.2:  Changed input from ".dat" to ".csv"
        Added example files in "/Check"

0.0.1:  Initial release.
'''

class Kite_Plot:
    '''
    Takes given distributions and creates a figure with corresponding kite plots.

    WARNING: At the moment the indices are used as x-values.
    '''
    def __init__(self) -> None:
        '''Initializes the kite_plot class.'''
        self.data_set = []               
        self.data_mean = []             
        self.data_std = []              
        self.offset = 0.0                
        self.lables = []                 

    def read_data(self, arr:list) -> None:
        '''Reads an array as a data set.'''
        self.data_set.append(arr)       
        self.data_mean.append(0)         
        self.data_std.append(0)          
        self.lables.append("x") 
        return         
    
    def calc_mean(self) -> None:
        '''Calculates the mean for every data set.'''
        for i, data in enumerate(self.data_set):
            mean = 0.0
            norm = 0.0
            for index,data_point in enumerate(data):
                mean += data_point*index
                norm += data_point
            mean /= norm
            self.data_mean[i] = mean
        return

    def calc_std(self) -> None:
        '''Calculates the standard deviation for every data set.'''
        self.calc_mean()
        for i, data in enumerate(self.data_set):
            std_dev = 0.0
            for index, data_point in enumerate(data):
                std_dev += ((self.data_mean[i] - index)**2)*data_point
            std_dev = std_dev**(0.5)
            self.data_std[i] = std_dev
        return

    def calc_offset(self) -> None:
        '''Calculates the offset necessary for the visualization as kite plot.'''
        max_value = 0.0
        for data in self.data_set:
            if max_value < max(data):
                max_value = max(data)
        self.offset = max_value*2.0
        return

    def calc_width(self,y_value:float,x_new:float,arr:list) -> list:
        '''Calculates the width of the kite plot at a specific y-value.'''
        idx = (np.abs(x_new - y_value)).argmin()
        return arr[idx]

    def write_plot(self,path:str="Kite_Plot.png") -> None:
        '''Plots the distributions as a figure with kite plots.'''
        plt.figure(figsize=(12,10))

        colors = ['darkblue',
                  'cornflowerblue',
                  'darkgreen',
                  'lightgreen',
                  'gold',
                  'yellow',
                  'red',
                  'lightcoral']

        self.calc_std()
        self.calc_offset()

        xticks = []

        for i, data in enumerate(self.data_set):
            color = colors[i % len(colors)]
            x = np.arange(len(data))
            x_new = np.linspace(x.min(),x.max(),300)
            spl = make_interp_spline(x,data,k=3)
            data_smooth = spl(x_new)

            offset_step = i*self.offset
            xticks.append(offset_step)

            plt.plot(data_smooth + offset_step, x_new, color, 
                     lw = 2, alpha = 0.6)
            plt.plot(-data_smooth + offset_step, x_new, color, 
                     lw = 2, alpha = 0.6)
            plt.fill_betweenx(x_new,-data_smooth+offset_step,data_smooth+offset_step,
                              color=color,alpha=0.2)
            
            mean_width = self.calc_width(self.data_mean[i],x_new,data_smooth)
            plt.plot([-mean_width+offset_step,mean_width+offset_step],
                     [self.data_mean[i],self.data_mean[i]], 
                     color="k",linestyle="--", linewidth=2)
            
            std_dev_plus = self.calc_width(self.data_mean[i]+self.data_std[i],
                                           x_new, data_smooth)
            std_dev_minus = self.calc_width(self.data_mean[i]-self.data_std[i],
                                           x_new, data_smooth)
            plt.plot([-std_dev_plus+offset_step,std_dev_plus+offset_step],
                     [self.data_mean[i]+self.data_std[i],self.data_mean[i]+self.data_std[i]],
                     color="k",linestyle=":",linewidth=1.5)
            plt.plot([-std_dev_minus+offset_step,std_dev_minus+offset_step],
                     [self.data_mean[i]-self.data_std[i],self.data_mean[i]-self.data_std[i]],
                     color="k",linestyle=":",linewidth=1.5)
            
        # Hier kann die Beschriftung des Grafen geändert werden
        plt.yticks(fontsize=16)
        plt.xticks(xticks,self.lables,rotation=45,ha='right',fontsize=16)

        plt.ylabel('Order Parameter',fontsize=16)

        plt.savefig(path)

        return

    def read_file(self,path:str) -> None:
        '''Reads the columns of the input file as data sets.'''
        columns = []
        with open(path,newline="") as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                print(row[1])
                for i,value in enumerate(row[1:]):
                    if (len(columns)) <= i:
                        columns.append([])
                    columns[i].append(float(value))
        for column in columns:
            self.read_data(column)
        csvfile.close()
        return
                
    def read_lables(self,arr:list) -> None:
        '''Takes an array as discription for the kites.'''
        for index,name in enumerate(arr):
            self.lables[index] = name
        return


def main(files:list,lables:list,path:str="Kite_Plot.png"):
    '''Plots the data of a given file and label list.'''
    kite = Kite_Plot()
    for file in files:
        kite.read_file(file)
    kite.read_lables(lables)
    kite.write_plot(path)



files = ["0/order.csv"]

lables = ["MC",
          "MD"]

main(files,lables)