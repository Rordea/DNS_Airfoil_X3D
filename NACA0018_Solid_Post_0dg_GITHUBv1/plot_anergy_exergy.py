# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 22:59:07 2024

@author: ricar
"""
import os
import matplotlib.pyplot as plt

from matplotlib import rcParams
import pandas as pd
import xcompact3d_toolbox as x3d


rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['font.size'] = 26


prm = x3d.Parameters(loadfile="./input.i3d")


def find_word(word_to_find):
    with open("./input.i3d", 'r') as file:
        for line in file:
            if word_to_find in line:
                words = line.split()
                if len(words) >= 3:
                    try:
                        third_word_float = float(words[2])
                        return third_word_float
                    except ValueError:
                        print("Third word is not a valid float.")
                else:
                    print("Line does not have enough words to extract the third word.")
                break
        else:
            print(f"Word '{word_to_find}' not found in any line.")
            return None






cv_x_min = find_word("!xld(1)")
cv_x_max = find_word("!xrd(1)")
cv_y_min = find_word("!yld(1)")
cv_y_max = find_word("!yud(1)")
cv_z_min = 0
cv_z_max = find_word("zlz ")

cv_x_max_0=cv_x_min+1

solid_airfoil_file = './results_data/Anergy_Exergy.dat'
porous_airfoil_file = '../NACA0018_Prs_Post_0dg/results_data/Anergy_Exergy.dat'

# Read the data, ensuring that the last column 'Cd total' is not split
df_solid_airfoil = pd.read_csv(solid_airfoil_file, delim_whitespace=True)
df_porous_airfoil = pd.read_csv(porous_airfoil_file, delim_whitespace=True)

# Extract coordinates and variables
coordinates_solid = df_solid_airfoil.iloc[:, 0]
variables_solid = df_solid_airfoil.iloc[:, 1:]

coordinates_porous = df_porous_airfoil.iloc[:, 0]
variables_porous = df_porous_airfoil.iloc[:, 1:]

# Plot 1: Coordinates vs. Variables from the solid airfoil file
plt.figure(figsize=(10, 6))

plt.plot(coordinates_solid, variables_solid.iloc[:, 0], label="$\epsilon_{a}$",color="blue")
plt.plot(coordinates_solid, variables_solid.iloc[:, 1], label="$\epsilon_{v}$",color="orange")
plt.plot(coordinates_solid, variables_solid.iloc[:, 2], label="$\epsilon_{p}$",color="green")
#plt.plot(coordinates_solid, variables_solid.iloc[:, 3], label="Etau")
#plt.plot(coordinates_solid, variables_solid.iloc[:, 4], label="Ek")
plt.plot(coordinates_solid, variables_solid.iloc[:, 5], label="Anergy ($\phi$)",color="red")
plt.plot(coordinates_solid, variables_solid.iloc[:, 6], label="Exergy ($\epsilon_{total}$)",color="purple")
plt.plot(coordinates_solid, variables_solid.iloc[:, 7], label="Cd total$_{ARF} (\phi+\epsilon_{total})$",color="brown")

plt.axvline(x=cv_x_min, color='red', linestyle=':', linewidth=1, alpha=0.9)
plt.axhline(y=0.1409, color='black', linestyle='--',label="Cd$_{RRF}$")
plt.axvline(x=cv_x_max_0, color='red', linestyle=':', linewidth=1, alpha=0.9, label="Airfoil limits")

plt.title('Energy Analysis \n Solid Airfoil')
plt.xlabel("x")
plt.ylabel('Cd')
#plt.legend(loc='best', ncol=2, frameon=False)
plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=4, frameon=False)
plt.tight_layout()
plt.minorticks_on()
plt.grid(which='both', color='gray', linestyle='-', linewidth=0.5, alpha=0.2)
plt.grid(True)
plt.show()

# Plot 2: Coordinates vs. Variables from the porous airfoil file


# Plot 2: Coordinates vs. Variables from the solid airfoil file
plt.figure(figsize=(10, 6))

plt.plot(coordinates_porous, variables_porous.iloc[:, 0], label="$\epsilon_{a}$",color="blue")
plt.plot(coordinates_porous, variables_porous.iloc[:, 1], label="$\epsilon_{v}$",color="orange")
plt.plot(coordinates_porous, variables_porous.iloc[:, 2], label="$\epsilon_{p}$",color="green")
#plt.plot(coordinates_solid, variables_solid.iloc[:, 3], label="Etau")
#plt.plot(coordinates_solid, variables_solid.iloc[:, 4], label="Ek")
plt.plot(coordinates_porous, variables_porous.iloc[:, 5], label="Anergy ($\phi$)",color="red")
plt.plot(coordinates_porous, variables_porous.iloc[:, 6], label="Exergy ($\epsilon_{total}$)",color="purple")
plt.plot(coordinates_porous, variables_porous.iloc[:, 7], label="Cd total$_{ARF} (\phi+\epsilon_{total})$",color="brown")

plt.axvline(x=cv_x_min, color='red', linestyle=':', linewidth=1, alpha=0.9)
plt.axhline(y=0.1409, color='black', linestyle='--',label="Cd$_{RRF}$")
plt.axvline(x=cv_x_max_0, color='red', linestyle=':', linewidth=1, alpha=0.9, label="Airfoil limits")

plt.title('Energy Analysis \n Porous Airfoil')
plt.xlabel("x")
plt.ylabel('Cd')
plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=4, frameon=False)
plt.tight_layout()
plt.minorticks_on()
plt.grid(which='both', color='gray', linestyle='-', linewidth=0.5, alpha=0.2)
plt.grid(True)
plt.show()

#Plot 3: Comparison of Variables between the solid and porous airfoil files
plt.figure(figsize=(10, 6))

plt.plot(coordinates_solid, variables_solid.iloc[:, 0], label="$\epsilon_{a}$",color="blue")
plt.plot(coordinates_solid, variables_solid.iloc[:, 1], label="$\epsilon_{v}$",color="orange")
plt.plot(coordinates_solid, variables_solid.iloc[:, 2], label="$\epsilon_{p}$",color="green")
#plt.plot(coordinates_solid, variables_solid.iloc[:, 3], label="Etau")
#plt.plot(coordinates_solid, variables_solid.iloc[:, 4], label="Ek")
plt.plot(coordinates_solid, variables_solid.iloc[:, 5], label="Anergy ($\phi$)",color="red")
plt.plot(coordinates_solid, variables_solid.iloc[:, 6], label="Exergy ($\epsilon_{total}$)",color="purple")
plt.plot(coordinates_solid, variables_solid.iloc[:, 7], label="Cd total$_{ARF} (\phi+\epsilon_{total})$",color="brown")






plt.plot(coordinates_porous, variables_porous.iloc[:, 0], label="$\epsilon_{a}$",color="blue",linestyle="--")
plt.plot(coordinates_porous, variables_porous.iloc[:, 1], label="$\epsilon_{v}$",color="orange",linestyle="--")
plt.plot(coordinates_porous, variables_porous.iloc[:, 2], label="$\epsilon_{p}$",color="green",linestyle="--")
#plt.plot(coordinates_solid, variables_solid.iloc[:, 3], label="Etau")
#plt.plot(coordinates_solid, variables_solid.iloc[:, 4], label="Ek")
plt.plot(coordinates_porous, variables_porous.iloc[:, 5], label="Anergy ($\phi$)",color="red",linestyle="--")
plt.plot(coordinates_porous, variables_porous.iloc[:, 6], label="Exergy ($\epsilon_{total}$)",color="purple",linestyle="--")
plt.plot(coordinates_porous, variables_porous.iloc[:, 7], label="Cd total$_{ARF} (\phi+\epsilon_{total})$",color="brown",linestyle="--")

plt.axvline(x=cv_x_min, color='red', linestyle=':', linewidth=1, alpha=0.9)
plt.axhline(y=0.1409, color='black', linestyle='--',label="Cd$_{RRF}$")
plt.axvline(x=cv_x_max_0, color='red', linestyle=':', linewidth=1, alpha=0.9, label="Airfoil limits")

plt.title('Energy Analysis \n Solid (-) vs Porous (--)')
plt.xlabel("x")
plt.ylabel('Cd')
plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=4, frameon=False)
plt.tight_layout()
plt.minorticks_on()
plt.grid(which='both', color='gray', linestyle='-', linewidth=0.5, alpha=0.2)
plt.grid(True)
plt.show()

