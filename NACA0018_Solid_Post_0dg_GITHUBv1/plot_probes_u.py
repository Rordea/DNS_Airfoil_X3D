# -*- coding: utf-8 -*-
"""
Created on Wed May 29 13:02:08 2024

@author: ricar
"""

# Filename: read_dat_file.py
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os

# Filename: read_and_plot_dat_file.py

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.size'] = 20

def read_dat_file(filename):
    data = {
        'column1': [],
        'column2': [],
        'column3': [],
        'column4': [],
        
    }

    with open(filename, 'r') as file:
        for line in file:
            # Split the line into three columns
            columns = line.split()
            
            # Ensure there are exactly three columns
            if len(columns) != 4:
                print(f"Skipping line: {line.strip()} (does not have exactly three columns)")
                continue

            # Append the values to the respective columns
            data['column1'].append(float(columns[0]))
            data['column2'].append(float(columns[1]))
            data['column3'].append(float(columns[2]))
            data['column4'].append(float(columns[3]))

    return data

def plot_columns(data,column1,column2,coeff,colors,file,markertype,linestyles):
    #plt.figure(figsize=(10, 6))
    plt.plot(column1, column2, label=f'{coeff}',color=colors,marker=markertype,markevery=20000,linestyle=linestyles)

    

# Example usage:
if __name__ == "__main__":
    
    probes=["probe0001","probe0002","probe0003","probe0004"]
    markertype=["v","x","+","o"]
    linestyles=["-",":","-.","--"]
    vels=["u","v","w"]
    for f in probes:
        
        
        file=f
        filename = f'./results_data/{file}'  # Replace with your .dat file name
        data = read_dat_file(filename)
        
        colors=["red","black","blue","purple","green"]
        
        coeff = "Ux"
        column1 = data['column1']
        column2 = data['column2']  
        plot_columns(data,column1,column2,vels[0],colors[0],file,markertype[0],linestyles[0])
        
        
        coeff = "Uy"
        column1 = data['column1']
        column2 = data['column3']
        plot_columns(data,column1,column2,vels[1],colors[1],file,markertype[1],linestyles[1])
        
        coeff = "Uz"
        column1 = data['column1']
        column2 = data['column4']
        plot_columns(data,column1,column2,vels[2],colors[2],file,markertype[2],linestyles[2])
    
    
        plt.xlabel('Time')
        plt.ylabel("u,v,w")
        plt.title(f'{f} \n velocities vs Time')
        plt.legend()
        #plt.xlim(0,6)
        #plt.ylim(0.0001,0.000125)
        plt.minorticks_on()
        plt.grid(which='both', color='gray', linestyle='-', linewidth=0.5, alpha=0.2)
        plt.grid(True)
        plt.legend(framealpha=0.0)
        plt.tight_layout()
        plt.savefig(f"./{f}_vel_vs_Time_{file}.png")
        plt.show()      

    os.system("pause")