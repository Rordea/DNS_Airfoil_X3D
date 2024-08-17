"""
       FILE: cylinder.py
     AUTHOR: Paul Bartholomew
DESCRIPTION: Post processes the cylinder case with reference data.
"""
import numpy as np
import matplotlib.pyplot as plt
from py4incompact3d.postprocess.postprocess import Postprocess
import xcompact3d_toolbox as x3d
from matplotlib.colors import LogNorm
import matplotlib.colors as mcolors
import os


def pause():
    input("Press Enter to continue...")



#;*****;*****;*****;*****;*****;*****;*****;*****;*****;***** P R E P A R I N G     .J S O N        F I L E  ;*****;*****;*****;*****;*****;*****;*****    
file_input3d = os.path.join(os.path.abspath(os.path.join(os.getcwd(), os.pardir)), 'input.i3d')

def deactivate_linux_vars(var_to_disable):
    
    file_path = file_input3d    ;     word_to_find = var_to_disable
    
    
    
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
handling_variables = ["ilist", "wrotation", "spinup_time",  "nstat","initstat", "xyzprobes", "iforces", "nvol", "xld",   "xrd",        "yld",    "yud",  "flag_all_digits",  "flag_extra_probes"]
[deactivate_linux_vars(f) for f in handling_variables]

prm = x3d.Parameters(loadfile=file_input3d)
nx_jsn = prm.nx        ; ny_jsn = prm.ny            ;  nz_jsn = prm.nz
lx_jsn = prm.xlx       ; ly_jsn = prm.yly           ; lz_jsn = prm.zlz
bcx_jsn = prm.nclx1    ; bcy_jsn = prm.ncly1        ; bcz_jsn=prm.nclz1
beta_jsn = prm.beta    ; strtch_jsn = prm.istret    ; timestep = prm.icheckpoint

def replace_line(file_path, variable_name, new_line):
    updated_lines = []
    variable_found = False

    with open(file_path, 'r') as file:
        for line in file:
            # Strip leading/trailing whitespace for comparison
            #stripped_line = line.strip()

            # Check if the line starts with the variable_name
            if line.startswith(variable_name):
                # Replace the entire line with the new_line
                updated_lines.append(new_line + '\n')
                variable_found = True
            else:
                # Keep the line unchanged
                updated_lines.append(line)

    if not variable_found:
        print(f"{variable_name} not found in the file.")
        return False

    # Write the updated content back to the file
    with open(file_path, 'w') as file:
        file.writelines(updated_lines)
    
    return True

file_path = os.path.join(os.getcwd(), 'input.json')
print("preparing input.json file")
variable_name = '      "Nx"' ; new_line = f'      "Nx": {nx_jsn},' ; replace_line(file_path, variable_name, new_line)
variable_name = '      "Ny"' ; new_line = f'      "Ny": {ny_jsn},' ; replace_line(file_path, variable_name, new_line)
variable_name = '      "Nz"' ; new_line = f'      "Nz": {nz_jsn},' ; replace_line(file_path, variable_name, new_line)
variable_name = '      "Lx"' ; new_line = f'      "Lx": {lx_jsn},' ; replace_line(file_path, variable_name, new_line)
variable_name = '      "Ly"' ; new_line = f'      "Ly": {ly_jsn},' ; replace_line(file_path, variable_name, new_line)
variable_name = '      "Lz"' ; new_line = f'      "Lz": {lz_jsn},' ; replace_line(file_path, variable_name, new_line)
variable_name = '      "BCx"' ; new_line = f'      "BCx": {bcx_jsn},' ; replace_line(file_path, variable_name, new_line)
variable_name = '      "BCy"' ; new_line = f'      "BCy": {bcy_jsn},' ; replace_line(file_path, variable_name, new_line)
variable_name = '      "BCz"' ; new_line = f'      "BCz": {bcz_jsn},' ; replace_line(file_path, variable_name, new_line)
variable_name = '          "beta"' ; new_line = f'          "beta": {beta_jsn},' ; replace_line(file_path, variable_name, new_line)
variable_name = '          "stretched"' ; new_line = f'          "stretched": {strtch_jsn}' ; replace_line(file_path, variable_name, new_line)
print("input.json file is ready")
pause()
#;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****





chord=150 # [mm]
width=25 # [mm]
height=27 # [mm]

prm = x3d.Parameters(loadfile="../input.i3d")
loc_perc_x = 23 
airfl_init_x = prm.xlx*(loc_perc_x/100)

AoA=0
min_boundx=0
airfl_init_x+(min_boundx/chord)

leading_edge_location = airfl_init_x
trailing_edge_location= airfl_init_x+1




#;*****;*****;*****;*****;*****;*****;*****;*****;*****;***** U S E R   I N P U T  ;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****  

plt.rcParams.update({    'font.size': 18,})

INPUT="input.json"     ;         last_timestep=str(timestep)
location_Slice_XZ = 0.5 # PROPORTION
location_Slice_YZ = 0.5 # PROPORTION
#Profile along Y
x_c = [0.2,0.4,0.6,0.8,1,1.2]  # METERS max. 4 points
location_x_0 = []
for f in x_c:
    print(f)
    print(f+leading_edge_location)
    location_x_0.append(f+leading_edge_location)
    
location_z_0 = 0.08   # METERS max. 4 points
print(location_x_0)
pause()
#;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****











#;*****;*****;*****;*****;*****;*****;*****;*****;*****;***** F U N C T I O N S  ;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****  

def read_mesh(timestep):
    postprocess = Postprocess(INPUT)
    mesh = postprocess.mesh        ;          t = timestep ###TO MODIFY TO MODIFY###
    postprocess.load(time=[t])
    return mesh
    

def read_variable(timestep,var=""):
    
    postprocess = Postprocess(INPUT)     ;    t = timestep ###TO MODIFY TO MODIFY###
    postprocess.load(time=[t])
    var_s = postprocess.fields[f"{var}"].data[t]     
    return var_s


def get_coordinates():
    x_coordinates=mesh.get_grid()[0] ; y_coordinates=mesh.get_grid()[1] ; z_coordinates=mesh.get_grid()[2]
    return x_coordinates , y_coordinates , z_coordinates



cool_to_warm_extended_colors = [
    (59 / 255, 76 / 255, 192 / 255),    (68 / 255, 90 / 255, 204 / 255),    (77 / 255, 104 / 255, 215 / 255),    (87 / 255, 117 / 255, 225 / 255),    (98 / 255, 130 / 255, 234 / 255),    (108 / 255, 142 / 255, 241 / 255),
    (119 / 255, 154 / 255, 247 / 255),    (130 / 255, 165 / 255, 251 / 255),    (141 / 255, 176 / 255, 254 / 255),    (152 / 255, 185 / 255, 255 / 255),    (163 / 255, 194 / 255, 255 / 255),    (174 / 255, 201 / 255, 253 / 255),
    (184 / 255, 208 / 255, 249 / 255),    (194 / 255, 213 / 255, 244 / 255),    (204 / 255, 217 / 255, 238 / 255),    (213 / 255, 219 / 255, 230 / 255),    (221 / 255, 221 / 255, 221 / 255),    (229 / 255, 216 / 255, 209 / 255),
    (236 / 255, 211 / 255, 197 / 255),    (241 / 255, 204 / 255, 185 / 255),    (245 / 255, 196 / 255, 173 / 255),    (247 / 255, 187 / 255, 160 / 255),    (247 / 255, 177 / 255, 148 / 255),    (247 / 255, 167 / 255, 135 / 255),
    (244 / 255, 157 / 255, 123 / 255),    (241 / 255, 147 / 255, 111 / 255),    (236 / 255, 136 / 255, 99 / 255),    (229 / 255, 125 / 255, 88 / 255),    (222 / 255, 114 / 255, 77 / 255),    (213 / 255, 103 / 255, 66 / 255),
    (203 / 255, 92 / 255, 56 / 255),    (192 / 255, 80 / 255, 46 / 255),    (180 / 255, 68 / 255, 37 / 255),    (167 / 255, 56 / 255, 29 / 255),    (152 / 255, 42 / 255, 22 / 255),    (136 / 255, 27 / 255, 16 / 255),
    (118 / 255, 8 / 255, 13 / 255),]




def cool_to_warm_extended():
    return mcolors.LinearSegmentedColormap.from_list(
        "cool_to_warm_extended", cool_to_warm_extended_colors, N=256)


def slice_xz(loc):
    
    if isinstance(loc, int) or isinstance(loc, float):
        print("Locating slice according to the proportion established")
        dim1 = (mesh.Lz/(1/loc))  /mesh.dz
        
    elif isinstance(loc, str):
        print("Locating slice according to the number of nodes established")
        dim1 = (loc)  
    slice_var = umean[:, :, int(dim1)]   
    return slice_var

    
def slice_yz(loc):
    
    if isinstance(loc, int) or isinstance(loc, float):
        print("Locating slice according to the proportion established")
        dim1 = (mesh.Lx/(1/loc))  /mesh.dx
        
    elif isinstance(loc, str):
        print("Locating slice according to the number of nodes established")
        dim1 = (loc)  

    slice_var = umean[int(dim1), :, :]
    return slice_var

def plot_slice_variable(a_coordinates,b_coordinates,slice_z0,name,slice_format="normal"):
    
    selected_cmap=cool_to_warm_extended()

    if slice_format=="normal":
        
        X, Y = np.meshgrid(a_coordinates, b_coordinates, indexing='ij')
        plt.contourf(X, Y, slice_z0, cmap=selected_cmap, levels = 256)
        plt.colorbar(label='Normalized Umean')
        plt.savefig(f"umean{name}.png",bbox_inches="tight")
        plt.show()
        
    elif slice_format=="inverted":
        Y, X = np.meshgrid(a_coordinates, b_coordinates, indexing='ij')
        plt.contourf(X, Y, slice_z0, cmap=selected_cmap, levels = 250)
        plt.colorbar(label='Normalized Umean')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.savefig(f"umean{name}.png",bbox_inches="tight")
        plt.show()

def plot_profile_of_variable_in_y(coordinates, location_x, location_z, variable, name,colors,linestyles,markers):
    
    dimX = location_x/mesh.dx
    dimZ = location_z/mesh.dz
    line_var = variable[int(dimX), :, int(dimZ)]
    
    with open('vel_profiles.dat', mode='w') as file:
        # Write the header for all columns
        location_x_0_str=[]
        for f in x_c:
            location_x_0_str.append(str(f))
        
        header = ['y-coords'] + location_x_0_str
        file.write('\t\t '.join(header) + '\n')
    
    
        for i in range(len(b_coordinates)):
            # Start the row with the current y-coordinate
            row = [b_coordinates[i]]
            
            # Append data for each profile
            for location_x in location_x_0:
                dimX = location_x / mesh.dx
                dimZ = location_z / mesh.dz
                line_var = variable[int(dimX), :, int(dimZ)]
                row.append(line_var[i])
            
            
            row_string = '\t\t'.join(map(str, row))
            file.write(row_string + '\n')
        
    print("velocity profiles have been written to vel_prfiles.dat")
    pause()
        # # Write each row of data
        # num_rows = len(b_coordinates)
        # for i in range(num_rows):
        #     # Create a row starting with the y-coordinate
        #     row = [b_coordinates[i]]
            
        #     # Append data for each column
        #     for name in location_x_0:
        #         row.append(umean[name][i])
    
        #     # Convert each element to a string and join with spaces
        #     row_string = ' '.join(map(str, row))
        #     file.write(row_string + '\n')
            
    plt.plot(line_var,coordinates , label=f'Umean at (x={location_x:.2f}, z={location_z:.2f})', color=colors,linestyle=linestyles,marker=markers,markevery=50)
    plt.xlabel(f'{name} ')
    plt.ylabel('Y coordinates of the Domain')
     
def plot_slice_variable_and_mesh(a_coordinates,b_coordinates,slice_z0,name,zoomXmin,zoomXmax,zoomYmin,zoomYmax):

    X, Y = np.meshgrid(a_coordinates, b_coordinates, indexing='ij')

    plt.contourf(X, Y, slice_z0, cmap=cool_to_warm_extended(), levels = 250)
    for i in range(len(b_coordinates)):
        plt.plot(a_coordinates, [b_coordinates[i]] * len(a_coordinates), color='black', linewidth=0.5)    
    
    # Plot vertical lines
    for i in range(len(a_coordinates)):
        plt.plot([a_coordinates[i]] * len(b_coordinates), b_coordinates, color='black', linewidth=0.5)
    
    plt.ylim(zoomYmin,zoomYmax)
    plt.xlim(zoomXmin,zoomXmax)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig("mesh.png",bbox_inches="tight")
    plt.show()


def find_word(word_to_find):
    with open(file_input3d, 'r') as file:
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



def volume_of_control_airfoil(mesh):
    
    dim1x = int(cv_x_min / mesh.dx)
    dim2x = int(cv_x_max / mesh.dx)
    dim1z = int(cv_z_min / mesh.dz)
    dim2z = int(cv_z_max / mesh.dz)
    dim1y = 232  
    dim2y = 816
    
    return dim1x, dim2x, dim1y, dim2y, dim1z, dim2z

def get_coordinates_for_vol(x0, x1, y0, y1, z0, z1, mesh):
    x_coordinates = mesh.get_grid()[0][x0:x1]
    y_coordinates = mesh.get_grid()[1][y0:y1]
    z_coordinates = mesh.get_grid()[2][z0:z1]
    return x_coordinates, y_coordinates, z_coordinates

#;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****














#;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****R E A D I N G    C O N T R O L    V O L U M E  ;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****;*****  

cv_x_min = find_word("!xld(1)")
cv_x_max = find_word("!xrd(1)")
cv_y_min = find_word("!yld(1)")
cv_y_max = find_word("!yud(1)")
cv_z_min = 0
cv_z_max = find_word("zlz ")

print("\n\n\n==================================================================")
print("Volume Control limits from input.i3d \n")
print(f"xmin = {cv_x_min}")
print(f"xmax = {cv_x_max}")
print(f"ymin = {cv_y_min}")
print(f"ymax = {cv_y_max}")
print(f"zmin = {cv_z_min}")
print(f"zmax = {cv_z_max}")
print("========================================================================")

pause()



    # Create the plot
mesh = read_mesh(last_timestep)
x0, x1, y0, y1, z0, z1 = volume_of_control_airfoil(mesh)
a_coordinates, b_coordinates, c_coordinates=get_coordinates_for_vol(x0, x1, y0, y1, z0, z1 , mesh)

# VARIABLE
umean = read_variable(last_timestep,var="umean")
vmean = read_variable(last_timestep,var="vmean")
#wmean = read_variable(last_timestep,"wmean")
pmean = read_variable(last_timestep,var="pmean")


mesh = read_mesh(last_timestep)





# ;****;****;****;****;****;****; I N L E T  # ;****;****;****;****;****;****;
node_numbers = np.arange(y0, y1-y0)  
#output_file = 'NP_CVinlet_data'
output_file = 'NP_CVoutlet_data'
with open(output_file, 'w') as file:
    file.write(f"{'node':>10}    {'x-coordinate':>14}      {'y-coordinate':>14}     {'dy-velocity-dy':>16} "
               f"   {'dy-velocity-dx':>16}    {'dx-velocity-dy':>16}   {'dx-velocity-dx':>16} "
               f"   {'y-velocity':>12}             {'x-velocity':>12}                {'pressure':>12}\n")
    for i in range(1,len(node_numbers)):
 
        file.write(f"        {node_numbers[i]}    {cv_x_min}    {b_coordinates[i]}    {dy_velocity_dy[x0][i][int(z1/2)]}      {dy_velocity_dx[x0][i][int(z1/2)]}             {dx_velocity_dy[x0][i][int(z1/2)]}            {dx_velocity_dx[x0][i][int(z1/2)]}        {vmean[x0][i][int(z1/2)]}           {umean[x0][i][int(z1/2)]}            {pmean[x0][i][int(z1/2)]} \n")

print(f"Data has been written to {output_file}")


#;****;****;****;****;****;****;;****;****;****;****;****;****;;****;****;****;****;****;****;;****;****;****;****;****;****;;****;****;****;****;****;****;






# ;****;****;****;****;****;****; O U T L E T  # ;****;****;****;****;****;****;
node_numbers = np.arange(y0, y1-y0)  
output_file = 'NP_CVinlet_data'
#output_file = 'NP_CVoutlet_data'
with open(output_file, 'w') as file:
    file.write(f"{'node':>10}    {'x-coordinate':>14}      {'y-coordinate':>14}     {'dy-velocity-dy':>16} "
               f"   {'dy-velocity-dx':>16}    {'dx-velocity-dy':>16}   {'dx-velocity-dx':>16} "
               f"   {'y-velocity':>12}             {'x-velocity':>12}                {'pressure':>12}\n")
    for i in range(1,len(node_numbers)):
 
        file.write(f"        {node_numbers[i]}    {cv_x_max}    {b_coordinates[i]}    {dy_velocity_dy[x1][i][int(z1/2)]}      {dy_velocity_dx[x1][i][int(z1/2)]}             {dx_velocity_dy[x1][i][int(z1/2)]}            {dx_velocity_dx[x1][i][int(z1/2)]}        {vmean[x1][i][int(z1/2)]}           {umean[x1][i][int(z1/2)]}            {pmean[x1][i][int(z1/2)]} \n")

print(f"Data has been written to {output_file}")
#;****;****;****;****;****;****;;****;****;****;****;****;****;;****;****;****;****;****;****;;****;****;****;****;****;****;;****;****;****;****;****;****;




# ;****;****;****;****;****;****; T O P  # ;****;****;****;****;****;****;
node_numbers = np.arange(x0, x1-x0)  
output_file = 'NP_CVtop_data'
with open(output_file, 'w') as file:
    file.write(f"{'node':>10}    {'x-coordinate':>14}      {'y-coordinate':>14}     {'dy-velocity-dy':>16} "
               f"   {'dy-velocity-dx':>16}    {'dx-velocity-dy':>16}   {'dx-velocity-dx':>16} "
               f"   {'y-velocity':>12}             {'x-velocity':>12}                {'pressure':>12}\n")
    for i in range(1,len(node_numbers)):
 
        file.write(f"        {node_numbers[i]}    {a_coordinates[i]}    {cv_y_max}    {dy_velocity_dy[i][y1][int(z1/2)]}      {dy_velocity_dx[i][y1][int(z1/2)]}             {dx_velocity_dy[i][y1][int(z1/2)]}            {dx_velocity_dx[i][y1][int(z1/2)]}        {vmean[i][y1][int(z1/2)]}           {umean[i][y1][int(z1/2)]}            {pmean[i][y1][int(z1/2)]} \n")

print(f"Data has been written to {output_file}")
#;****;****;****;****;****;****;;****;****;****;****;****;****;;****;****;****;****;****;****;;****;****;****;****;****;****;;****;****;****;****;****;****;







# ;****;****;****;****;****;****; B O T T O M  # ;****;****;****;****;****;****;
node_numbers = np.arange(x0, x1-x0)  
output_file = 'NP_CVbottom_data'
with open(output_file, 'w') as file:
    file.write(f"{'node':>10}    {'x-coordinate':>14}      {'y-coordinate':>14}     {'dy-velocity-dy':>16} "
               f"   {'dy-velocity-dx':>16}    {'dx-velocity-dy':>16}   {'dx-velocity-dx':>16} "
               f"   {'y-velocity':>12}             {'x-velocity':>12}                {'pressure':>12}\n")
    for i in range(1,len(node_numbers)):
 
        file.write(f"        {node_numbers[i]}    {a_coordinates[i]}    {cv_y_min}    {dy_velocity_dy[i][y0][int(z1/2)]}      {dy_velocity_dx[i][y0][int(z1/2)]}             {dx_velocity_dy[i][y0][int(z1/2)]}            {dx_velocity_dx[i][y0][int(z1/2)]}        {vmean[i][y0][int(z1/2)]}           {umean[i][y0][int(z1/2)]}            {pmean[i][y0][int(z1/2)]} \n")

print(f"Data has been written to {output_file}")
#;****;****;****;****;****;****;;****;****;****;****;****;****;;****;****;****;****;****;****;;****;****;****;****;****;****;;****;****;****;****;****;****;






# Horizontal Plane
a_coordinates , b_coordinates, c_coordinates = get_coordinates()

# slice_z0=slice_xz(location_Slice_XZ)
# name="_XZ_Horz_slc_"+str(location_Slice_XZ)
# plot_slice_variable(a_coordinates, b_coordinates,slice_z0,name,slice_format="normal")

# zoomXmin= 2.28 ; zoomXmax= 3.3 ; zoomYmin=4.89 ; zoomYmax = 5.10
# plot_slice_variable_and_mesh(a_coordinates,b_coordinates,slice_z0,name,zoomXmin,zoomXmax,zoomYmin,zoomYmax,)

# zoomXmin= 0 ; zoomXmax= 10 ; zoomYmin=0 ; zoomYmax = 10
# plot_slice_variable_and_mesh(a_coordinates,b_coordinates,slice_z0,name,zoomXmin,zoomXmax,zoomYmin,zoomYmax,)

# # Vertical Plane

# slice_z0=slice_yz(location_Slice_YZ)
# name="_YZ_Vrtcl_slc_"+str(location_Slice_YZ)
# plot_slice_variable(b_coordinates, c_coordinates,slice_z0,name,slice_format="inverted")


colors=["blue","black","red","purple","orange","green","cyan","orange","maroon","gray","lightgreen","pink"]
linestyles=["-","--",":","-.","-.",":","--","-"]
markers=["x","+","v","o","x","+","x","v","o","x"]

for f in range(0,len(location_x_0)):
    plot_profile_of_variable_in_y(b_coordinates, location_x_0[f], location_z_0, umean,name="umean",colors=colors[f],linestyles=linestyles[f],markers=markers[f])
    
    



    

    
plt.title('Umean')
plt.legend()
plt.minorticks_on()
plt.grid(which='both', linestyle='--', linewidth=0.5)

plt.savefig("Umean_line_example.png", bbox_inches="tight")
plt.show()


















#    file_path = "./umean.txt"

# Open the file in write mode ('w' mode)
#    with open(file_path, 'w') as file:
#        file.write(str(umean))




input("Press Enter to continue...")
#umean = avg_over_axis(mesh, umean, 2)