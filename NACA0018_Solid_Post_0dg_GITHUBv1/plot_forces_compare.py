# -*- coding: utf-8 -*-
"""
Created on Wed May 29 13:02:08 2024

@author: ricar
"""

# Filename: read_dat_file.py
import matplotlib.pyplot as plt
import os

# Filename: read_and_plot_dat_file.py

import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np




rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.size'] = 20




filename_baseline = './results_data/forces.dat' 

filename_to_compare = os.path.join(os.path.abspath(os.path.join(os.getcwd(), os.pardir)), 'NACA0018_Prs_Post_0dg/results_data/forces.dat' )



Time_start_avg = 7
color_solid="darkred"
color_porous="green"
avg_every = 50




def read_dat_file(filename):
    data = {
        'column1': [],
        'column2': [],
        'column3': []
    }

    with open(filename, 'r') as file:
        for line in file:
            # Split the line into three columns
            columns = line.split()
            
            # Ensure there are exactly three columns
            if len(columns) != 3:
                print(f"Skipping line: {line.strip()} (does not have exactly three columns)")
                continue

            # Append the values to the respective columns
            data['column1'].append(float(columns[0]))
            data['column2'].append(float(columns[1]))
            data['column3'].append(float(columns[2]))

    return data

def plot_columns_solid(data,column1,column2,coeff,name):
    #plt.figure(figsize=(10, 6))
    plt.plot(column1, column2, '-', label=f'{name}',color=color_solid)

    
def plot_columns_porous(data,column1,column2,coeff,name):
    #plt.figure(figsize=(10, 6))
    plt.plot(column1, column2, '-', label=f'{name}',color=color_porous)



def read_data(coeff,data): 
    
    global Time, coeff_aero
    
    if coeff.startswith("Cd"):
        Time = data['column1']
        coeff_aero = data['column2']
    
    
    elif coeff.startswith("CL"):
        Time = data['column1']
        coeff_aero = data['column3']
    

def find_closest_index(time_array, target_value):
    time_np=np.array(time_array)
    differences = np.abs(time_np - target_value)
    closest_index = np.argmin(differences)
    return closest_index

# Example usage:
if __name__ == "__main__":
    
    coefficients = ["Cd","CL"]
    
    for f in range(0,len(coefficients)):
    
        coeff = coefficients[f]
    
        data_baseline = read_dat_file(filename_baseline)   
        read_data(coeff, data_baseline)
        
        index = find_closest_index(Time, Time_start_avg)
        
        index_fin = find_closest_index(Time, Time[-1])
        
        avg_time=[]
        Time_mod=[]
        std_dev = []
        convergence_abs=[]
        
        for f in range(index,index_fin,avg_every):
            print(f"Averaging... {index_fin-f} timesteps to go")
            # print(coeff_aero[index])
            # print(coeff_aero[index:f+1])
            # print(f-index+2)
            
            avg = (coeff_aero[index]+np.sum(coeff_aero[index:f+1]))/(f-index+2)
            avg_time.append(avg)
            std_dev.append(np.std(avg_time,ddof=1))
            convergence_abs.append(coeff_aero[f+1]-coeff_aero[f])
            Time_mod.append(Time[f])
            
        
        data_to_compare = read_dat_file(filename_to_compare)
        read_data(coeff, data_to_compare)
        
        avg_time1=[]
        Time_mod1=[]
        std_dev1 = []
        convergence_abs1=[]
        for f in range(index,index_fin,avg_every):
            print(f"Averaging... {index_fin-f} timesteps to go")
            # print(coeff_aero[index])
            # print(coeff_aero[index:f+1])
            # print(f-index+2)
            
            avg = (coeff_aero[index]+np.sum(coeff_aero[index:f+1]))/(f-index+2)
            avg_time1.append(avg)
            std_dev1.append(np.std(avg_time,ddof=1))
            convergence_abs1.append(coeff_aero[f+1]-coeff_aero[f])
            Time_mod1.append(Time[f])
            
            
            
        print(f"\n{coeff}_avg solid : {avg_time[-1]}")
        print(f"\n{coeff}_avg porous : {avg_time1[-1]}")
        
        print(f"{coeff}_inst porous : {coeff_aero[-1]}")
        
        data_baseline = read_dat_file(filename_baseline)   
        read_data(coeff, data_baseline)
        print(f"{coeff}_inst solid : {coeff_aero[-1]}")
        
        print(f"\nStandard deviation of the moving average of the {coeff} solid        : {std_dev[-1]}")
        print(f"Average of the absolute errors of {coeff} solid                        : {np.mean(convergence_abs)}") 
        
        print(f"\nStandard deviation of the moving average of the {coeff} porous       : {std_dev1[-1]}")
        print(f"Average of the absolute errors of {coeff} porous                       : {np.mean(convergence_abs1)}")
        
        os.system("pause")
        
        data_baseline = read_dat_file(filename_baseline)   
        read_data(coeff, data_baseline)         
        plt.figure(figsize=(12, 8))
        plot_columns_solid(data_baseline,Time,coeff_aero,coeff,"Solid 0 deg")
        
        data_to_compare = read_dat_file(filename_to_compare)
        read_data(coeff, data_to_compare)
        
        coeff_aero1=coeff_aero
        
        plot_columns_porous(data_to_compare,Time,coeff_aero1,coeff,"Porous 0 deg")
        
        plt.plot(Time_mod,avg_time,label="Solid Mvg Avg",linestyle="--",color=color_solid,marker="x",markevery=avg_every*4)
        plt.plot(Time_mod1,avg_time1,label="Porous Mvg Avg",linestyle="--",color=color_porous,marker="v",markevery=avg_every*4)
        
        plt.axvline(x=Time_start_avg, color='black', linestyle='--', linewidth=1, alpha=1,label="Start Avg")
        plt.xlim(Time_start_avg-1,Time[-1]+0.5)
        plt.title(f'{coeff} vs Time',fontsize=40)
        plt.xlabel('Time',fontsize=40)
        plt.ylabel(f'{coeff}',fontsize=40)
        plt.minorticks_on()
        plt.tick_params(axis='both', which='both', width=2, length=2,labelsize=28)
        plt.gca().spines['top'].set_linewidth(1.5)
        plt.gca().spines['bottom'].set_linewidth(1.5)
        plt.gca().spines['left'].set_linewidth(1.5)
        plt.gca().spines['right'].set_linewidth(1.5)
 
        plt.grid(which='both', color='gray', linestyle='-', linewidth=1, alpha=0.2)
        
        adjusting_limits={coeff_aero[-1],coeff_aero1[-1],avg_time[-1],avg_time1[-1]}
        min_lim=min(  adjusting_limits  )
        
        max_lim=max(  adjusting_limits  )
        
        plt.ylim(min_lim-(20*std_dev[-1]),max_lim+(20*std_dev[-1]))
        #plt.legend(framealpha=0.0)
        legend = plt.legend( loc='best', ncol=2,fontsize=25,frameon=False)

        plt.tight_layout()

        plt.savefig(f"./{coeff}_vs_Time_comparison0_4dg.png")
        plt.show()
    
os.system("pause")