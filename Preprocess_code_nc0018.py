# Import necessary libraries
import xcompact3d_toolbox as x3d
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import math
import os
import shutil 
import zipfile
import re
import sys
import time





######################## INPUT PARAMETERS NEEDED FOR THIS CODE TO RUN #################################
def input():
    
    global naca,file_input3d,chord,width,height,max_boundy,min_boundx,min_boundy,AoA, \
    CAD,probes_in_holes,pressure_taps,points_line_probe,vol_contrl_threshold,\
    low_vol_contrl_threshold,Final_sim_time,output_images_time
    
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
    
    AoA= 0 # deg
    
    CAD = f'C:/Users/ricar/Documents/IRP_2/NACA_{naca}_CAD/wing_nc{naca}_{AoA}dg.stl'
    
    probes_in_holes = "no"
    pressure_taps = "no"
    
    points_line_probe = 5
    
    
    
    vol_contrl_threshold = 0.065
    low_vol_contrl_threshold = 0.015
    
    Final_sim_time = 10 # [ seconds ]
    output_images_time = 0.5 # [ seconds ]

# nraf each refinement layer takes 9[  minutes] aprox. with Beta = 1, 15 min. with Beta = 0.5


input()

# Pressure taps input

airfoil_file_text_upper = 'NACA0018_upper.txt'
airfoil_file_text_lower = 'NACA0018_lower.txt'
pressure_taps_offst_to_curve = 0.0065
#pressure_taps_offst_to_curve = 0.01
number_of_pressure_taps = 10
 


########################################################################################################

# I N P U T      O F     P R O B E S   I N S I D E      T H E        H O L E S
    
                           #InletX  #InletY  #InletZ
probes_line_coords_inlet = [(0.6962,      0.516,   0.0833), #1st Layer
                                                                       
                            (0.6962,      0.4999,   0.0833), #2nd Layer
                            
                            (0.6962,      0.4821,   0.0833), #3rd Layer
                            
                            (0.7207,      0.466,   0.0833), #4th Layer
                            
                            (0.7505,      0.449,   0.0833), #5th Layer
                            
                            (0.8012,      0.4326,   0.0833), #6th Layer
                            
                            (0.9316,      0.415,   0.0833),]#7th Layer





probes_line_coords_outlet = [(1.6043,      0.516,   0.0833), #1st Layer
                                                                       
                            (1.6699,      0.4999,   0.0833), #2nd Layer
                            
                            (1.5719,      0.4821,   0.0833), #3rd Layer
                            
                            (1.4783,      0.466,   0.0833), #4th Layer
                            
                            (1.3882,      0.449,   0.0833), #5th Layer
                            
                            (1.2684,      0.4329,   0.0833), #6th Layer
                            
                            (1.0584,      0.4153,   0.0833),]#7th Layer







###################################################################################################













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
    print(point)
    print(rotated_point)
    
    
    
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







origin_airfoil=[airfl_init_x+(chord/chord),airfl_init_y+((height/chord)/2),loc_init_z+((width/chord)/2)]




axis = origin_airfoil[0], origin_airfoil[1], origin_airfoil[2]





#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,   PROBES IN HOLES (IDENTIFICATION AND ROTATION )  **,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,




if probes_in_holes=="yes":

    for f  in range(0,len(probes_line_coords_inlet)):  
        for g in range(0,2,1):
            
            if g == 0:
                deltax = probes_line_coords_outlet[f][g] - probes_line_coords_inlet[f][g]            
                deltay = probes_line_coords_outlet[f][g+1] - probes_line_coords_inlet[f][g+1]
                z_coord = probes_line_coords_outlet[f][g+2]
            else:
                print("")
        
            m = deltay/deltax
            
            
        

     
        dxx=(probes_line_coords_outlet[f][0]-probes_line_coords_inlet[f][0])/(points_line_probe+1)
        for n in range(1,points_line_probe+1):
             xx = probes_line_coords_inlet[f][0]+(dxx*n)
             p = m*(xx) + probes_line_coords_outlet[f][1]
             
             
             point = xx , p, origin_airfoil[2]
             axis = origin_airfoil[0], origin_airfoil[1], origin_airfoil[2]   # Rotate around the z-axis
             angle = np.radians(AoA)  
    
             rotated_point = rotate_point(point, axis, angle)
             x_rot, y_rot ,z_rot = rotated_point
             
             # print(f"{rotated_point}")
             # print(f"{x_rot}, {y_rot} , {z_rot}")
             
            
             
             # print(f"point {n} coordinates : ({xx},{p},{z_coord})")
             
             probes.append(rotated_point)


else:
    print("\nno probes in holes\n")

#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,     













start_time = time.perf_counter()

check_data_folder_to_exist("./data", directory_path_zip)


print(f"{int(xlx_txt)}")
os.system("pause")

check_data_folder_to_exist(f"./run_{int(xlx_txt)}_{prm.nx}_{int(yly_txt)}_{int(prm.ny)}_{int(zlz_txt)}_{int(prm.nz)}_{prm.beta}.zip", directory_path_zip)







#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,   DEFINING CONTROL VOLUME TO MEASURE FORCES      )  **,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,

#

vol_contr_x_left = (airfl_init_x-0.01)
vol_contr_x_right = (airfl_init_x+1.03)
vol_contr_y_bottom = (loc_centr_y-(height/chord)/2)-low_vol_contrl_threshold
vol_contr_y_top = (loc_centr_y+(height/chord)/2)+vol_contrl_threshold+max_boundy


# WRITING CONTROL VOLUME

word_to_find = '!xld(1)'    # Replace with the word you want to find
new_phrase = f'!xld(1) = {vol_contr_x_left} !X left for volume control'  # Replace with the new phrase
replace_phrases_in_input_file(word_to_find, new_phrase)

word_to_find = '!xrd(1)'    # Replace with the word you want to find
new_phrase = f'!xrd(1) = {vol_contr_x_right} !X left for volume control'  # Replace with the new phrase
replace_phrases_in_input_file(word_to_find, new_phrase)

word_to_find = '!yld(1)'    # Replace with the word you want to find
new_phrase = f'!yld(1) = {vol_contr_y_bottom} !Y bottom for volume control'  # Replace with the new phrase
replace_phrases_in_input_file(word_to_find, new_phrase)

word_to_find = '!yud(1)'    # Replace with the word you want to find
new_phrase = f'!yud(1) = {vol_contr_y_top} !Y top for volume control'  # Replace with the new phrase
replace_phrases_in_input_file(word_to_find, new_phrase)


#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,















#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,  WRITING PROBES  **,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,


probes_text=[]
probes_text.append("! Probes in Domain \n")
for f in range(0,4,1):
    for j in range(0,3,1):
        
        probe_coord = probes[f][j]
        
        probe_text = f'!xyzprobes({j+1},{f+1}) = {probe_coord}    ! position of probe in domain {f+1}\n'  # Replace with the new phrase
        probes_text.append(probe_text)
probes_text.append("\n\n")
# WRITING PROBE LOCATIONS IN input.i3d 
if probes_in_holes=="yes":
   
    probes_text.append("! Probes in Holes \n")
    for f in range(4,(len(probes)),1):
        for j in range(0,3,1):
            
            probe_coord = probes[f][j]
            
            probe_text = f'!xyzprobes({j+1},{f+1}) = {probe_coord}    ! position of probe in holes {f+1}\n'  # Replace with the new phrase
            probes_text.append(probe_text)
    probes_text.append("\n\n")
    
elif pressure_taps=="yes":
   
    probes_text.append("! Pressure taps \n")
    for f in range(4,(len(probes)),1):
        for j in range(0,3,1):
            
            probe_coord = probes[f][j]
            
            probe_text = f'!xyzprobes({j+1},{f+1}) = {probe_coord}    ! position of pressure tap {f+1}\n'  # Replace with the new phrase
            probes_text.append(probe_text)
    probes_text.append("\n\n")
    
    
word1 = '&ProbesParam'
word2 = '&ScalarParam'
delete_and_add_probes_between_words(file_input3d, word1, word2, probes_text)






#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,





























#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**, READ STL FILE   **,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,

origin_point = airfl_init_y
adj_leg = axis[0]-(airfl_init_x+0.3) # adjust maximum thickness position for better results in case there is other airfoil.
opsste_leg = adj_leg*(math.tan(np.radians(AoA)))

#xer=0.01
#yer=0.0167

#xer=0.001428*AoA # 7 deg
#yer=0.002385*AoA # 7 deg

#move the point with 0.001 of error

#xer=0.001428*AoA # 4 deg
#yer=0.001385*AoA # 4 deg


xer=0*AoA # 0 deg
yer=0*AoA # 0 deg




x3d.param["mytype"] = np.float64
epsi = x3d.init_epsi(prm, dask=True)

file_path = os.path.join(os.path.abspath(os.path.join(os.getcwd(), os.pardir)), CAD) #finds the file in the parent folder
for key in epsi.keys():
    #print(key)                                                                                        #this number adjusts the final y coord of the geometry with the probes.
    epsi[key] = epsi[key].geo.from_stl(file_path, origin=dict(x=airfl_init_x+(min_boundx/chord),y=axis[1]+(min_boundy/chord),z=airfl_init_z),scale=1.0/chord, 
                                       rotate=dict(axis=[0,0,1],theta=math.radians(0)), 
                                       user_tol=1e-5)
    
    print(epsi[key])
    
dataset = x3d.gene_epsi_3D(epsi, prm)

#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,






















#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**, PLOT HORIZONTAL PLANE, VERTICAL PLANE WITH CONTROL VOLUME AND DOMAIN PROBES, HOLES PROBES, PRESSURE TAPS   **,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,

# PLOT HORIZONTAL PLANE
epsi["epsi"].sel(z=loc_init_z+((width/chord)/2), method="nearest").plot(x="x",cbar_kwargs={'orientation': 'horizontal'})
end_time = time.perf_counter()
colors="white"


# PLOT CONTROL VOLUME

plt.plot([vol_contr_x_left, vol_contr_x_right], [vol_contr_y_bottom, vol_contr_y_bottom], linestyle='--',color=colors)
plt.plot([vol_contr_x_left, vol_contr_x_right], [vol_contr_y_bottom, vol_contr_y_bottom], linestyle='--',color=colors)

plt.plot([vol_contr_x_left, vol_contr_x_right], [vol_contr_y_top, vol_contr_y_top], linestyle='--',color=colors)
plt.plot([vol_contr_x_left, vol_contr_x_right], [vol_contr_y_top, vol_contr_y_top], linestyle='--',color=colors)

plt.plot([vol_contr_x_left, vol_contr_x_left], [vol_contr_y_top, vol_contr_y_bottom], linestyle='--',color=colors)
plt.plot([vol_contr_x_left, vol_contr_x_left], [vol_contr_y_top, vol_contr_y_bottom], linestyle='--',color=colors)

plt.plot([vol_contr_x_right, vol_contr_x_right], [vol_contr_y_top, vol_contr_y_bottom], linestyle='--',color=colors)
plt.plot([vol_contr_x_right, vol_contr_x_right], [vol_contr_y_top, vol_contr_y_bottom], linestyle='--',color=colors)


# PLOT PROBES

for f in range(0,4,1):
    
    plt.plot(probes[f][0],probes[f][1] , marker="x",color="red")
    plt.annotate(
    f'Probe {f+1}',      
    xy=(probes[f][0], probes[f][1]), 
    xytext=(probes[f][0] + 0.015, probes[f][1] + 0.005), 
    fontsize=8,           # Font size
    color='white',           # Text color
    )

plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("latest_mesh_representation_z0plane_test.png")
plt.show()


# PLOT PROBES IN HOLES OR PRESSURE TAPS OR BOTH

if probes_in_holes=="yes" or pressure_taps=="yes":
    
    for f in range(4,(len(probes)),1):
        
        plt.plot(probes[f][0],probes[f][1] , marker="x",color=colors)
        plt.annotate(
        f'Probe {f+1}',      
        xy=(probes[f][0], probes[f][1]), 
        xytext=(probes[f][0] + 0.015, probes[f][1] + 0.005), 
        fontsize=8,           # Font size
        color='black',           # Text color

        )
    
epsi["epsi"].sel(z=loc_init_z+((width/chord)/2), method="nearest").plot(x="x",cbar_kwargs={'orientation': 'horizontal'})
plt.gca().set_aspect('equal', adjustable='box')
plt.plot(axis[0],axis[1],color="white",marker="+",markersize=5)
plt.ylim(vol_contr_y_bottom+0.01,vol_contr_y_top+0.1)
plt.xlim(vol_contr_x_left+0.01,vol_contr_x_right+0.01)
plt.savefig("latest_mesh_representation_z0plane_zoom.png")
plt.show()
    
    
    
    
    
    
    
epsi["epsi"].sel(x=2.6, method="nearest").plot(y="y",cbar_kwargs={'orientation': 'vertical'})
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("latest_mesh_representation_x0plane_test.png")

plt.show()    

epsi["epsi"].sel(x=2.6, method="nearest").plot(y="y",cbar_kwargs={'orientation': 'vertical'})
plt.gca().set_aspect('equal', adjustable='box')
plt.ylim(vol_contr_y_bottom+0.01,vol_contr_y_top-0.37)
plt.xlim(0,0.1666666666667)
plt.savefig("latest_mesh_representation_x0plane_test_zoom.png")

plt.show()    
 
#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,



















#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**, WRITING DATA FOLDER   **,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,



if prm.iibm >= 2:
    prm.nobjmax = dataset.obj.size
ds  = x3d.init_dataset(prm)




# Set boundary conditions
for key in 'bxx1 bxy1 bxz1'.split():

    ds[key] *= 0.0
    if key=='bxx1':
        ds[key] += 1.0



# Set initial conditions
for key in 'ux uy uz '.split():

    ds[key] *= 0.0
    if key=='ux':
        ds[key] += 1.0

# for key in 'vol_frc'.split():
#     print(ds[key].attrs['name'])
#     ds[key] *= 0.0
#     # if key=='ux':
#     #     ds[key] += 1.0

# Write the fields
prm.dataset.write(ds)



#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,











# ACTIVATING VARIABLES AGAIN TO USE THEM IN XCOMPACT3D SIMULATION

for f in handling_variables:
    "\n\nActivating variables.... \n\n"
    var_to_activate=f
    activate_linux_vars(var_to_activate)




word_to_find = 'nprobes'    # Replace with the word you want to find
new_phrase = f'nprobes={len(probes)}'  # Replace with the new phrase
replace_phrases_in_input_file(word_to_find, new_phrase)

word_to_find = 'ilast'   
new_phrase = f'ilast = {timesteps}'  
replace_phrases_in_input_file(word_to_find, new_phrase) 

word_to_find = 'icheckpoint'   
new_phrase = f'icheckpoint = {int((timesteps-1)/2)}'  
replace_phrases_in_input_file(word_to_find, new_phrase) 

word_to_find = 'iprocessing'   
new_phrase = f'iprocessing = {int((timesteps-1))}'  
replace_phrases_in_input_file(word_to_find, new_phrase) 

word_to_find = 'initstat'   
new_phrase = f'initstat = {int((timesteps-1)/2)}'  
replace_phrases_in_input_file(word_to_find, new_phrase) 

word_to_find = 'ioutput'   
new_phrase = f'ioutput = {int((output_images_time)/prm.dt)}'  
replace_phrases_in_input_file(word_to_find, new_phrase) 










total_time= end_time-start_time

print(f"\n\nElapsed time: {total_time:.6f} seconds\n\n")

print(f"\n\nElapsed time: {((total_time)/60):.6f} minutes\n\n")

print("Zipping file.... ")

folder_to_zip = files_to_zip
zip_file_name = directory_path_zip

zip_files(folder_to_zip, zip_file_name)


print("\n\ndata.zip file and input.i3d are ready to upload in HPC\n\n")

print("\n\nTo submit to HPC:\n"
      "Upload data.zip and input.i3d to your directory and run the script run.sh.\n"
      "Make sure you have the mpi_2022.sub file and the run.sh in the directory.")





os.system("pause")