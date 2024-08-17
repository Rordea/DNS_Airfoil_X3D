# Import necessary libraries
import xcompact3d_toolbox as x3d
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pressure_taps import press_taps
import math
import os
import shutil 
import zipfile
import re
import sys
import time





######################## INPUT PARAMETERS NEEDED FOR THIS CODE TO RUN #################################

naca = "0018"

file_input3d = "input.i3d"


chord=150 # [mm]
width=25 # [mm]
height=27 # [mm]

max_boundy = 49.0299/chord #at 17 deg. #max lim in the volume of control. Dont move


min_boundy = -13.5045

min_boundx = 0
#max_boundx = 1000.17
#max_boundy = 72.5441

AoA= 4 # deg

CAD = f'C:/Users/ricar/Documents/IRP_2/NACA_{naca}_CAD/wing_nc{naca}_{AoA}dg.stl'

probes_in_holes = "no"
pressure_taps = "yes"

points_line_probe = 5



vol_contrl_threshold = 0.065
low_vol_contrl_threshold = 0.015

Final_sim_time = 10 # [ seconds ]
output_images_time = 0.5 # [ seconds ]

# nraf each refinement layer takes 9[  minutes] aprox. with Beta = 1, 15 min. with Beta = 0.5




# Pressure taps input

airfoil_file_text_upper = 'NACA0018_upper.txt'
airfoil_file_text_lower = 'NACA0018_lower.txt'
pressure_taps_offst_to_curve = 0
#pressure_taps_offst_to_curve = 0.01
number_of_pressure_taps = 120
 














# START OF DEFINED FUNCTIONS

def zip_files(files_to_zip, zip_file_name):
    base_path = os.path.commonpath(files_to_zip)  # Common base path for relative paths
    with zipfile.ZipFile(zip_file_name, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for file in files_to_zip:
            if os.path.isfile(file):
                zipf.write(file, os.path.relpath(file, base_path))
            elif os.path.isdir(file):
                for root, dirs, files in os.walk(file):
                    for file_in_dir in files:
                        file_path = os.path.join(root, file_in_dir)
                        zipf.write(file_path, os.path.relpath(file_path, base_path))
                
                
                
def replace_phrases_in_input_file(word_to_find,new_phrase):
    
    file_path = file_input3d    # Replace with your file path
    word_to_find = word_to_find    # Replace with the word you want to find
    new_phrase = new_phrase
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Open the file in write mode to overwrite it with the modified content
    with open(file_path, 'w') as file:
        for f in range(0, len(lines)):
            # Skip empty lines
            if not lines[f].strip():
                file.write(lines[f])
                continue

            words_space = lines[f].split()
           
            
            # Check if the first word matches word_to_find
            if words_space[0] == word_to_find:
                print(f"Original line: {lines[f].strip()}")
                # Replace the entire line with new_phrase
                file.write(new_phrase + '\n')
                print(f"Modified lines[f]: {new_phrase}")
            else:
                # Write the original line if no replacement is needed
                file.write(lines[f])


def delete_and_add_probes_between_words(filename, word1, word2, new_lines):
 
    with open(filename, 'r') as file:
        lines = file.readlines()

    line_num1 = -1
    line_num2 = -1

    # Find the first occurrence of word1 and word2
    for i, line in enumerate(lines):
        if word1 in line and line_num1 == -1:
            line_num1 = i
        if word2 in line and line_num1 != -1:  # start looking for word2 only after finding word1
            line_num2 = i
            break

    # Check if both words were found
    if line_num1 == -1 or line_num2 == -1:
        print("Both words were not found in the file.")
        return

    # Delete lines between line_num1 and line_num2
    
    new_file_content = lines[:line_num1+3] + new_lines + lines[line_num2-4:]

    # Write the modified lines back to the file
    with open(filename, 'w') as file:
        file.writelines(new_file_content)

    print(f"Lines between {word1} and {word2} have been deleted and replaced with new lines.")




                
                
def deactivate_linux_vars(var_to_disable):
    
    file_path = file_input3d    # Replace with your file path
    word_to_find = var_to_disable    # Replace with the word you want to find
    
    
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    
    with open(file_path, 'w') as file:
        for line in lines:
            #print(line)
            
            if word_to_find in line:
                if not line.startswith("!"):
                    
    
                    modified_line = "!" + line
    
                    file.write(modified_line)
                elif line.startswith("!"):
    
                    file.write(line)
            else:
                file.write(line)
                
                
def activate_linux_vars(var_to_disable):
    
    file_path = file_input3d    # Replace with your file path
    word_to_find = var_to_disable    # Replace with the word you want to find
    
    
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    
    with open(file_path, 'w') as file:
        for line in lines:
            #print(line)
            
            if word_to_find in line:
                if not line.startswith("!"):
                    print(f"variable '{word_to_find}' already activated, no need to change anything") 
                    file.write(line)
                    
                elif line.startswith("!"):
                    print(f"Original line: {line.strip()}")
                    modified_line = line[1:]
                    print(f"the modified_line is {modified_line}")
                    
                    print(f"Activating variable {word_to_find} in line: {modified_line.strip()}")
                    file.write(modified_line)
            else:
                file.write(line)



def check_data_folder_to_exist(directory_path,directory_path_zip):


    directory_path = directory_path
    directory_path_zip = directory_path_zip
    
    # Attempt to remove the directory
    try:
        shutil.rmtree("./data")
        print(f"Directory '{directory_path}' removed successfully.")
    except FileNotFoundError:
        print(f"Directory '{directory_path}' not found. Continuing...")
    
    
    
    try:
        os.remove("./data.zip")
        print(f"Directory '{directory_path_zip}' removed successfully.")
    except FileNotFoundError:
        print(f"Directory '{directory_path_zip}' not found. Continuing...")
        
        
        


def rotate_point(point, axis, angle):
    """
    Rotate a point around a given axis by a given angle.
    
    :param point: A tuple or list of (x, y, z) coordinates.
    :param axis: A tuple or list of (x, y, z) axis coordinates to rotate around.
    :param angle: The angle in radians to rotate the point.
    :return: The rotated point as a NumPy array.
    """
    # Ensure the axis is a unit vector
    axis = np.array(axis)
    axis = axis / np.linalg.norm(axis)
    
    # Create the rotation matrix using the axis and angle
    cos_angle = np.cos(-angle)
    sin_angle = np.sin(-angle)
    one_minus_cos = 1.0 - cos_angle

    ux, uy, uz = axis
    rotation_matrix = np.array([
        [cos_angle + ux*ux*one_minus_cos, ux*uy*one_minus_cos - uz*sin_angle, ux*uz*one_minus_cos + uy*sin_angle],
        [uy*ux*one_minus_cos + uz*sin_angle, cos_angle + uy*uy*one_minus_cos, uy*uz*one_minus_cos - ux*sin_angle],
        [uz*ux*one_minus_cos - uy*sin_angle, uz*uy*one_minus_cos + ux*sin_angle, cos_angle + uz*uz*one_minus_cos]
    ])

    # Rotate the point
    rotated_point = np.dot(rotation_matrix, np.array(point))
    return rotated_point






     
# END OF DEFINED FUNCTIONS


def rotate_press_taps(point, axis_point, angle):
    """
    Rotate a point around the z-axis located at a specific point.
    
    :param point: A tuple or list of (x, y, z) coordinates.
    :param axis_point: The point around which to rotate (x, y, z).
    :param angle: The angle in radians to rotate the point.
    :return: The rotated point as a NumPy array.
    """
    # Translate the point to the origin
    translated_point = np.array(point) - np.array(axis_point)
    
    
    # Create the rotation matrix for the z-axis
    cos_angle = np.cos(-angle)
    sin_angle = np.sin(-angle)
    rotation_matrix = np.array([
        [cos_angle, -sin_angle, 0],
        [sin_angle, cos_angle, 0],
        [0, 0, 1]
    ])
    
    # Rotate the point
    rotated_translated_point = np.dot(rotation_matrix, translated_point)
    
    # Translate the point back
    rotated_point = rotated_translated_point + np.array(axis_point)
    #print(point)
    #print(rotated_point)
    
    
    
    return rotated_point
      
      
      
      
      
      
      
      
      
        
#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,   START  OF   THE  CODE     #**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,     




# Deactivate variables that xcompact3D-toolbox doesn't recognize but are needed in 
handling_variables = ["ilist", "wrotation", "spinup_time",  "nstat","initstat", "xyzprobes", "iforces", "nvol",
                        "xld",   "xrd",        "yld",    "yud",  "flag_all_digits",  "flag_extra_probes"]

for f in handling_variables:
    deactivate_linux_vars(f)


prm = x3d.Parameters(loadfile=file_input3d)

xlx_txt=prm.xlx
yly_txt=prm.yly
zlz_txt=prm.zlz


files_to_zip = [
    './input.i3d',
    './data',
    "./latest_mesh_representation_z0plane_zoom.png"
]

num_oder=prm.itimescheme
nraf=prm.nraf
#directory_path_zip = f"./run_x{xlx_txt:.0f}_{prm.nx}_y{yly_txt:.0f}_{prm.ny:.0f}_z{zlz_txt:.2f}_{prm.nz:.0f}_B{prm.beta:.2f}_Tmeord{num_oder}_nraf{nraf}.zip"
directory_path_zip = f"C:/Users/ricar/Documents/IRP_2/NACA_{naca}_CAD/{AoA}deg_run_x{xlx_txt:.0f}_{prm.nx}_y{yly_txt:.0f}_{prm.ny:.0f}_z{zlz_txt:.2f}_{prm.nz:.0f}_B{prm.beta:.2f}_Tmeord{num_oder}_nc{naca}.zip"


mesh_refined=prm.get_mesh()



x_mesh = mesh_refined["x"]
y_mesh = mesh_refined["y"]
z_mesh = mesh_refined["z"]

dx=[]
dy=[]
dz=[]

for f in range(1,len(x_mesh)):
    dxi=x_mesh[f]-x_mesh[f-1]
    dx.append(dxi)

for f in range(1,len(y_mesh)):
    dxi=y_mesh[f]-y_mesh[f-1]
    dy.append(dxi)
    
for f in range(1,len(z_mesh)):
    dxi=z_mesh[f]-z_mesh[f-1]
    dz.append(dxi)

min_dx=min(dx)
print(f"\n\nmin. edge length in x :  {min_dx}")
std_dx=np.std(dx)
print(f"standard deviation of the dx list : {std_dx} \n\n")


min_dy=min(dy)
print(f"\n\nmin. edge length in y :  {min_dy}")
std_dy=np.std(dy)
print(f"standard deviation of the dy list : {std_dy} \n\n")


min_dz=min(dz)
print(f"\n\nmin. edge length in z :  {min_dz}")
std_dz=np.std(dz)
print(f"standard deviation of the dz list : {std_dz} \n\n")


os.system("pause")




timesteps = int(Final_sim_time/prm.dt)

loc_perc_x = 23   # location of the starting point of the airfoil in percent of the X length of the domain.
loc_centr_y = prm.yly/2 #location of the center of the airfoil in y direction
loc_init_z = 0 # threshold of the airfoil to the domain in Z direction for both sides in CHORDS.




airfl_init_x = prm.xlx*(loc_perc_x/100)


airfl_init_y = loc_centr_y-((height/chord)/2)
airfl_init_z = loc_init_z

  
           # x    # y       # z
probes = [((0.1*prm.xlx),   (0.5*prm.yly),   ((width/chord)/2) ), # Probe 1
           
          (((airfl_init_x+(chord/chord))+(0.5*airfl_init_x)),   (0.5*prm.yly),   ((width/chord)/2)),    # Probe 2
          
          ((airfl_init_x+(0.25*chord/chord)),  (0.9*prm.yly),   ((width/chord)/2)),   # Probe 3
          
          (((airfl_init_x+(0.25*chord/chord))),  (0.1*prm.yly),   ((width/chord)/2)),]  # Probe 4 



#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,     













axis_point = (min_boundx+airfl_init_x+1, loc_centr_y, (width/chord)/2)  
angle = math.radians(AoA)
 

#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,   PRE VISUALISATION OF PRESSURE TAPS  **,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,


if pressure_taps=="yes":
     # Pressure Taps Upper
    press_taps_coords_upper = press_taps(filename=airfoil_file_text_upper, degrees=range(1, 17), offset_distance = pressure_taps_offst_to_curve ,    num_offset_points = number_of_pressure_taps   ,initial_point_x=airfl_init_x,initial_point_y=loc_centr_y,
                                   initial_point_z=((width/chord)/2))

    
        
    rotated_press_taps_coords_upper = []
    
    for point in press_taps_coords_upper['offset_points']:
        
        rotated_point = rotate_press_taps(point, axis_point, angle)
        rotated_press_taps_coords_upper.append(rotated_point)

    for i, point in enumerate(rotated_press_taps_coords_upper):
        print(f"Pressure tap {i+1} coordinates: ({point[0]:.2f}, {point[1]:.2f}, {point[2]:.2f})")
        #probes.append( (point[0], point[1], point[2]) )

    
        with open(f'Cp_coords_upper_{AoA}dg.dat', mode='w') as file:
            # Write the header for all columns
            
            header = ['x-coords'] + ['y-coords'] + ['z-coords']
            file.write('\t\t '.join(header) + '\n')
    
    
            # Write the data rows
            for i,point in enumerate(rotated_press_taps_coords_upper):
                row = [point[0],point[1],point[2]]
                
                file.write('\t\t'.join(map(str, row)) + '\n')

    x_fit = np.linspace(min(press_taps_coords_upper['offset_points'][:, 0]), max(press_taps_coords_upper['offset_points'][:, 0]), 100)
    y_fit_uppr = press_taps_coords_upper['best_polynomial'](x_fit)



    # Pressure Taps Lower
    press_taps_coords_lower = press_taps(filename=airfoil_file_text_lower, degrees=range(1, 17), offset_distance = -pressure_taps_offst_to_curve ,    num_offset_points = number_of_pressure_taps   ,initial_point_x=airfl_init_x,initial_point_y=loc_centr_y,
                                   initial_point_z=((width/chord)/2))
    
    
    
    
    
    
    
    
    
    rotated_press_taps_coords_lower = []
    
    for point in press_taps_coords_lower['offset_points']:
        rotated_point = rotate_press_taps(point, axis_point, angle)
        rotated_press_taps_coords_lower.append(rotated_point)
        
        
        
    # Plotting example (optional)
    
    for i, point in enumerate(rotated_press_taps_coords_lower):
       # print(f"Pressure tap {i+1} coordinates: ({point[0]:.2f}, {point[1]:.2f}, {point[2]:.2f})")
        
        
        
        with open(f'Cp_coords_lower_{AoA}dg.dat', mode='w') as file:
            # Write the header for all columns
            
            header = ['x-coords'] + ['y-coords'] + ['z-coords']
            file.write('\t\t '.join(header) + '\n')
    
    
            # Write the data rows
            for i,point in enumerate(rotated_press_taps_coords_lower):
                row = [point[0],point[1],point[2]]
                file.write('\t\t'.join(map(str, row)) + '\n')
        
        #probes.append( (point[0], point[1], point[2]) )
    #probes.append(press_taps_coords['offset_points'])




    x_fit = np.linspace(min(press_taps_coords_lower['offset_points'][:, 0]), max(press_taps_coords_lower['offset_points'][:, 0]), 100)
    y_fit_lwr = press_taps_coords_lower['best_polynomial'](x_fit)
    
    rotated_press_taps_coords_upper_array = np.array(rotated_press_taps_coords_upper)
    rotated_press_taps_coords_lower_array = np.array(rotated_press_taps_coords_lower)
    
    # print(rotated_press_taps_coords_upper_array.shape)
    # os.system("pause")
    
    plt.scatter(rotated_press_taps_coords_upper_array[:, 0], rotated_press_taps_coords_upper_array[:, 1], color='green', marker="+" ,label='Offset Points')
    plt.scatter(rotated_press_taps_coords_lower_array[:, 0], rotated_press_taps_coords_lower_array[:, 1], color='green', marker="+" ,label='Offset Points')
    
    plt.plot(x_fit, y_fit_lwr, label=f'Best Degree {press_taps_coords_upper["best_degree"]}, MSE: {press_taps_coords_upper["min_error"]:.2f}', color='blue')
    plt.plot(x_fit, y_fit_uppr, label=f'Best Degree {press_taps_coords_lower["best_degree"]}, MSE: {press_taps_coords_lower["min_error"]:.2f}', color='blue')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Offset Points and Best Polynomial Fit')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(f"./press_taps_{AoA}.png")
    plt.show()
            
elif pressure_taps=="no":
    print("no pressure taps included")
else:
    print("\n ERROR!!! Check if the pressure taps are considered or not ( only valid options are  yes or no cas sensitive )\n")
    os.system("pause") ; 
#     sys.exit()




#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,     




