"""
       FILE: cylinder.py
     AUTHOR: Paul Bartholomew
DESCRIPTION: Post processes the cylinder case with reference data.
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from py4incompact3d.postprocess.postprocess import Postprocess
import xcompact3d_toolbox as x3d
from matplotlib.colors import LogNorm
import matplotlib.colors as mcolors
from stl import mesh
from scipy.interpolate import interp2d
import math
from scipy.interpolate import RegularGridInterpolator
import gc

gc.disable()

plt.rcParams.update({
    'font.size': 10
})


#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,   U S E R    I N P U T  **,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,

INPUT="input.json"
input_file_path ="../input.i3d"

prm = x3d.Parameters(loadfile=input_file_path)

last_timestep="499998"

check_mesh_interpolation = "no"

#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,















#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,   DEFINED FUNCTIONS  **,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,

def pause():
    input("Press Enter to continue...")

def find_word(word_to_find):
    third_word_float = None
    with open(input_file_path, 'r') as file:
        for line in file:
            if word_to_find in line:
                words = line.split()
                if len(words) >= 3:
                    try:
                        third_word_float = float(words[2])
                    except ValueError:
                        print("Third word is not a valid float.")
                else:
                    print("Line does not have enough words to extract the third word.")
                break
        else:
            print("Word 'find_word' not found in any line.")
    return third_word_float

        
    with open(input_file_path, 'r') as file:
        # Read each line in the file
        for line in file:
            # Find the word you're interested in (for example, 'find_word')
            if word_to_find in line:
                # Split the line into words
                words = line.split()
                
                # Ensure the line has at least 3 words
                if len(words) >= 3:
                    # Get the third word (index 2) and convert to float
                    try:
                        third_word_float = float(words[2])
                       
                    except ValueError:
                        print("Third word is not a valid float.")
                else:
                    print("Line does not have enough words to extract the third word.")
                break  # Exit the loop after finding the word (if you only need one occurrence)
        else:
            print("Word 'find_word' not found in any line.")
            
        return third_word_float

def read_mesh(timestep):
    
    postprocess = Postprocess(INPUT)
    mesh = postprocess.mesh
    t = timestep ###TO MODIFY TO MODIFY###
    postprocess.load(time=[t])
    yp=mesh.stretched
    yp = mesh.get_grid()
    uumean = postprocess.fields["uumean"].data[t]
    
    return mesh
    
def read_variable(timestep,variable="umean"):
    
    postprocess = Postprocess(INPUT)
    t = timestep ###TO MODIFY TO MODIFY###
    postprocess.load(time=[t])
    umean = postprocess.fields[variable].data[t] 
    
    return umean

def get_coordinates_for_vol(x0,x1,y0,y1,z0,z1):
    x_coordinates=mesh.get_grid()[0][x0:x1]
    y_coordinates=mesh.get_grid()[1][y0:y1]
    z_coordinates=mesh.get_grid()[2][z0:z1] 

    return x_coordinates , y_coordinates , z_coordinates

def get_coordinates():
    x_coordinates=mesh.get_grid()[0]
    y_coordinates=mesh.get_grid()[1]
    z_coordinates=mesh.get_grid()[2]
    return x_coordinates , y_coordinates , z_coordinates

def cool_to_warm_extended():
    cool_to_warm_extended_colors = [
        (59 / 255, 76 / 255, 192 / 255),    (68 / 255, 90 / 255, 204 / 255),    (77 / 255, 104 / 255, 215 / 255),    (87 / 255, 117 / 255, 225 / 255),    (98 / 255, 130 / 255, 234 / 255),    (108 / 255, 142 / 255, 241 / 255),    (119 / 255, 154 / 255, 247 / 255),
        (130 / 255, 165 / 255, 251 / 255),    (141 / 255, 176 / 255, 254 / 255),    (152 / 255, 185 / 255, 255 / 255),    (163 / 255, 194 / 255, 255 / 255),    (174 / 255, 201 / 255, 253 / 255),    (184 / 255, 208 / 255, 249 / 255),    (194 / 255, 213 / 255, 244 / 255),
        (204 / 255, 217 / 255, 238 / 255),    (213 / 255, 219 / 255, 230 / 255),    (221 / 255, 221 / 255, 221 / 255),    (229 / 255, 216 / 255, 209 / 255),    (236 / 255, 211 / 255, 197 / 255),    (241 / 255, 204 / 255, 185 / 255),    (245 / 255, 196 / 255, 173 / 255),
        (247 / 255, 187 / 255, 160 / 255),    (247 / 255, 177 / 255, 148 / 255),    (247 / 255, 167 / 255, 135 / 255),    (244 / 255, 157 / 255, 123 / 255),    (241 / 255, 147 / 255, 111 / 255),    (236 / 255, 136 / 255, 99 / 255),    (229 / 255, 125 / 255, 88 / 255),
        (222 / 255, 114 / 255, 77 / 255),    (213 / 255, 103 / 255, 66 / 255),    (203 / 255, 92 / 255, 56 / 255),    (192 / 255, 80 / 255, 46 / 255),    (180 / 255, 68 / 255, 37 / 255),    (167 / 255, 56 / 255, 29 / 255),    (152 / 255, 42 / 255, 22 / 255),
        (136 / 255, 27 / 255, 16 / 255),    (118 / 255, 8 / 255, 13 / 255),
    ]
    return mcolors.LinearSegmentedColormap.from_list(
        "cool_to_warm_extended", cool_to_warm_extended_colors, N=256)       
        



def volume_of_control_airfoil(umean,vmean,wmean,pmean,uumean,uvmean,uwmean,vvmean,vwmean,wwmean):
    
    y = mesh.get_grid()[1][:] 
    dim1x = int(cv_x_min/mesh.dx)                ;     dim2x = int(cv_x_max/mesh.dx  )
    dim1z = int(cv_z_min/mesh.dz)                ;     dim2z = int((cv_z_max)/mesh.dz)    
    dim1y = (np.abs(y - cv_y_min)).argmin()      ;     dim2y = (np.abs(y - cv_y_max)).argmin()

    umean_control_vol = umean[dim1x:dim2x, dim1y:dim2y, dim1z:dim2z]
    vmean_control_vol = vmean[dim1x:dim2x, dim1y:dim2y, dim1z:dim2z]
    wmean_control_vol = wmean[dim1x:dim2x, dim1y:dim2y, dim1z:dim2z]
    pmean_control_vol = pmean[dim1x:dim2x, dim1y:dim2y, dim1z:dim2z]    
    
    uumean_control_vol = uumean[dim1x:dim2x, dim1y:dim2y, dim1z:dim2z]
    uvmean_control_vol = uvmean[dim1x:dim2x, dim1y:dim2y, dim1z:dim2z]
    uwmean_control_vol = uwmean[dim1x:dim2x, dim1y:dim2y, dim1z:dim2z]
    vvmean_control_vol = vvmean[dim1x:dim2x, dim1y:dim2y, dim1z:dim2z]
    vwmean_control_vol = vwmean[dim1x:dim2x, dim1y:dim2y, dim1z:dim2z]
    wwmean_control_vol = wwmean[dim1x:dim2x, dim1y:dim2y, dim1z:dim2z]

    return umean_control_vol, vmean_control_vol, wmean_control_vol , pmean_control_vol, uumean_control_vol ,uvmean_control_vol ,uwmean_control_vol ,vvmean_control_vol,vwmean_control_vol , wwmean_control_vol,dim1x ,dim2x , dim1y, dim2y , dim1z, dim2z, 

def get_centered_face_values(var, a_coords, b_coords, c_coords):

    
    # Ensure input coordinates are sorted and within the expected range
    if not (np.all(np.diff(a_coords) > 0) and np.all(np.diff(b_coords) > 0) and np.all(np.diff(c_coords) > 0)):
        raise ValueError("Input coordinates must be strictly increasing.")
    
    
    x_mid = (a_coords[:-1] + a_coords[1:]) / 2
    y_mid = (b_coords[:-1] + b_coords[1:]) / 2
    z_mid = (c_coords[:-1] + c_coords[1:]) / 2
    

    centered_values = np.zeros((len(a_coords), len(z_mid), len(y_mid)))
    
    for x_idx, x_val in enumerate(a_coords):
        var_slice = var[x_idx, :, :]
        interpolator = interp2d(c_coords, b_coords, var_slice, kind='linear')
        
        Z_mid, Y_mid = np.meshgrid(z_mid, y_mid, indexing='ij')
        interpolated_values = interpolator(z_mid, y_mid)
        centered_values[x_idx, :, :] = interpolated_values.T

    return centered_values, x_mid, y_mid, z_mid 


def get_centered_cell_values(var, a_coords, b_coords, c_coords):
    
    if not (np.all(np.diff(a_coords) > 0) and np.all(np.diff(b_coords) > 0) and np.all(np.diff(c_coords) > 0)):
        raise ValueError("Input coordinates must be strictly increasing.")
    
    x_mid = (a_coords[:-1] + a_coords[1:]) / 2
    y_mid = (b_coords[:-1] + b_coords[1:]) / 2
    z_mid = (c_coords[:-1] + c_coords[1:]) / 2
    
    X_mid, Y_mid, Z_mid = np.meshgrid(x_mid, y_mid, z_mid, indexing='ij')
    
    interpolator = RegularGridInterpolator((a_coords, b_coords, c_coords), var, method='linear')
    
    midpoints = np.array([X_mid.flatten(), Y_mid.flatten(), Z_mid.flatten()]).T
    
    centered_values = interpolator(midpoints).reshape(len(x_mid), len(y_mid), len(z_mid))
    
    return centered_values, x_mid, y_mid, z_mid



def get_areas(a_coords,b_coords):

    dz = np.diff(a_coords) 
    dy = np.diff(b_coords) 

    areas = np.outer(dy, dz)
    
    return areas,dy,dz

def get_volumes(a_coords,b_coords,c_coords):
    
    dx = np.diff(a_coords)
    dy = np.diff(b_coords) 
    dz = np.diff(c_coords) 
    
    dx_grid, dy_grid, dz_grid = np.meshgrid(dx, dy, dz, indexing='ij')
    volumes = dx_grid*dy_grid*dz_grid
    
    return volumes,dx,dy,dz



def plot_mesh(i,j,h="x",v="y",color="blue"):
    
    i, j = np.meshgrid(i, j)
    x_coords = i.flatten()
    y_coords = j.flatten()
    plt.scatter(x_coords, y_coords, c=color, marker='o',s=1)
    plt.savefig(f"./mesh_{h}{v}_cv.png")
    plt.show()
    

def plot_mesh_centered_values(i,j,k,L,h="x",v="y",color="blue",color_fc="red"):
    
    i, j = np.meshgrid(i, j)
    k, L = np.meshgrid(k, L)
    x_coords = i.flatten()
    y_coords = j.flatten()
    x_coords_fc = k.flatten()
    y_coords_fc = L.flatten()
    plt.scatter(x_coords, y_coords, c=color, marker='o',s=1)
    plt.scatter(x_coords_fc, y_coords_fc, c=color_fc, marker='x',s=1)
    #plt.savefig(f"./mesh_{h}{v}_cv.png")
    plt.show()


def integration_yz(var):
    yz_areas_na=yz_areas[np.newaxis,:,:]
    mult_areas = var * yz_areas_na
    
    
    
    sum_along_y=np.sum(mult_areas[:,:,:],axis=1)
    sum_along_z=np.sum(sum_along_y[:,:],axis=1)
    total_integral_zy=sum_along_z ;
    return total_integral_zy



def integration_xz(var):
    mult_areas = var[:-1,:,:] * xz_areas[:,np.newaxis,:]
    sum_along_z=np.sum(mult_areas[:,:,:],axis=2)
    sum_along_x=np.sum(sum_along_z[:,:],axis=0)
    total_integral_xz=sum_along_x ;
    return total_integral_xz


def compute_gradients_vel(velu,velv,velw,dx_vol,dy_vol,dz_vol):
    
    dudxyz = np.gradient(velu,dx_vol,dy_vol,dz_vol)
    dvdxyz = np.gradient(velv,dx_vol,dy_vol,dz_vol)
    dwdxyz = np.gradient(velw,dx_vol,dy_vol,dz_vol)
    
    
    return dudxyz, dvdxyz, dwdxyz  



def integration_vol(var,vols):
    mult_areas = var * vols
    total_vol_integral=np.sum(mult_areas)
    return total_vol_integral

#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,












 #**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,   READ FIELD VARIABLES  **,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,

umean = read_variable(last_timestep,variable="umean")
vmean = read_variable(last_timestep,variable="vmean")
wmean = read_variable(last_timestep,variable="wmean")
pmean = read_variable(last_timestep,variable="pmean")
uumean = read_variable(last_timestep,variable="uumean") 
uvmean = read_variable(last_timestep,variable="uvmean") 
uwmean = read_variable(last_timestep,variable="uwmean") 
vvmean = read_variable(last_timestep,variable="vvmean") 
vwmean = read_variable(last_timestep,variable="vwmean") 
wwmean = read_variable(last_timestep,variable="wwmean") 
mesh = read_mesh(last_timestep)


cv_x_min=find_word("xld(1)")              ;             cv_x_max_1=find_word("xrd(1)")+5  ; cv_x_max_0=find_word("xrd(1)")
   
nx_vc = int( (cv_x_max_1-cv_x_min)/prm.dx)


Ea_plot=[]
Ev_plot = []
Ep_plot = []
Etau_plot = []
Ek_plot = []
Anergy_plot = []
Exergy_plot=[]
Cd_tot_plot=[]
x_coords_plot=[]


# with open("energy_data.dat", 'w') as file:
# # Write the header
#     file.write("Ea_plot\t\t\tEv_plot\t\t\tEp_plot\t\t\tEtau_plot\t\t\tEk_plot\t\t\tAnergy_plot\t\t\tExergy_plot\t\t\tCd_tot_plot\t\t\tx_coords_plot\n")

#**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,
for f in range(nx_vc-3,0,-1):
    
    cv_x_min=find_word("xld(1)")              ;             cv_x_max=cv_x_max_1-((prm.dx)*f)
    cv_y_min = find_word("yld(1)")-4
    cv_y_max = find_word("yud(1)")+4
    cv_z_min=0                                ;             cv_z_max=find_word("zlz ")
    
    
    # print("\n\n\n==================================================================")
    # print("Volume Control limits from input.i3d ")
    # print("====================================================================================================================================")
    # print(f"xmin = {cv_x_min}")
    # print(f"xmax = {cv_x_max}")
    
    # print(f"ymin = {cv_y_min}")
    # print(f"ymax = {cv_y_max}")
    
    # print(f"zmin = {cv_z_min}")
    # print(f"zmax = {cv_z_max}")
    # print("==========================================================================================================================================")
        
    
    
    
    
    
    
    
    
    #**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,   COMPUTE GRADIENTS  **,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,**,
    
    u,  v , w , pp, uu,uv,uw,vv,vw,ww, x0 , x1 , y0, y1,z0,z1 = volume_of_control_airfoil(umean,vmean,wmean,pmean,uumean,uvmean,uwmean,vvmean,vwmean,wwmean)
    x , y , z = get_coordinates_for_vol(x0,x1,y0,y1,z0,z1)
    
    
    u_fc, x_mid, y_mid, z_mid=get_centered_face_values(u, x, y, z)
    v_fc, x_mid, y_mid, z_mid=get_centered_face_values(v, x, y, z)
    w_fc, x_mid, y_mid, z_mid=get_centered_face_values(w, x, y, z)
    p_fc, x_mid, y_mid, z_mid=get_centered_face_values(pp, x, y, z)
    
    # print("====================================================================================================================================")
    # print("Checking nodes in faces")
    # print("====================================================================================================================================")
    # print(f"\nx,y,z shapes of arrays with node values:   {x.shape},{y.shape}.{z.shape} \n")
    # print(f"x,y,z shapes of arrays with face-centered values:   {x_mid.shape},{y_mid.shape}.{z_mid.shape}")
    # print("====================================================================================================================================")
    nodes_mesh_color = "blue"
    face_centered_mesh_color = "red"
    
    if check_mesh_interpolation == "yes":
        
        plot_mesh(x, y,h="x",v="y",color=nodes_mesh_color)
        plot_mesh(z, y,h="z",v="y",color=nodes_mesh_color)
        plot_mesh_centered_values(x, y, x_mid, y_mid,h="x",v="y",color=nodes_mesh_color,color_fc=face_centered_mesh_color)
        plot_mesh_centered_values(z, y, z_mid, y_mid,h="z",v="y",color=nodes_mesh_color,color_fc=face_centered_mesh_color)
    
    elif check_mesh_interpolation =="no":
        
        print("Skipping check of interpolated values ")
    
    else: 
        print("\n\nCheck Mesh Option not recognized!!!!!! \n\n")
        print("Program Terminated\n")
        sys.exit()
    # print("\n\n\n====================================================================================================================================")
    # print("Checking centered magnitude values in faces")
    # print("====================================================================================================================================")
    # print(f"\n u,v,w shapes of arrays with node values            :     u={u.shape},  v={v.shape},  w={w.shape} \n")
    # print(f" u,v,w shapes of arrays with face-centered values   :   ufc={u_fc.shape},vfc={v_fc.shape},wfc={w_fc.shape}")
    
    #print(f"\nu values in first four nodes                                     {u[0,0,0]}      {u[1,0,0]}") 
    #print(f"u values in first three faces                                           {u_fc[0,0,0]:.6f}    {u_fc[1,0,0]:.6f}    {u_fc[2,0,0]:.6f}") 
    #print(f"verification of u values with an average in first three faces           {((u[0,0,0]+u[1,0,0]+u[0,1,0]+u[1,1,0])/4):.6f}    {((u[1,0,0]+u[2,0,0])/2):.6f}    {((u[2,0,0]+u[3,0,0])/2):.6f}") 
    #print("====================================================================================================================================")
    u_fc_reshape = u_fc.transpose(0, 2, 1)
    v_fc_reshape = v_fc.transpose(0, 2, 1)
    w_fc_reshape = w_fc.transpose(0, 2, 1)
    p_fc_reshape = p_fc.transpose(0, 2, 1)
    
    
    
    u_cc, x_mid_vol, y_mid_vol, z_mid_vol=get_centered_cell_values(u, x, y, z)
    v_cc, x_mid_vol, y_mid_vol, z_mid_vol=get_centered_cell_values(v, x, y, z)
    w_cc, x_mid_vol, y_mid_vol, z_mid_vol=get_centered_cell_values(w, x, y, z)
    p_cc, x_mid_vol, y_mid_vol, z_mid_vol=get_centered_cell_values(pp, x, y, z)
    
    volumes,dx,dy,dz=get_volumes(x, y, z)
    
    dx_vol=np.diff(x_mid_vol)
    dy_vol=np.diff(y_mid_vol)
    dz_vol=np.diff(z_mid_vol)
    
    # print("\n\n\n====================================================================================================================================")
    # print("Checking Centered Cell values procedure")
    # print("=======================================================================================================================================")
    # print(f"volumes array shape : {volumes.shape}")
    # # print(f"Last cell in volumes array {volumes[137,578,10]}")
    # # print(f"Last cell in volumes array {dx[137]*dy[578]*dz[10]}")
    
    # print(f"Cells U value array shape : {u_cc.shape}")
    # print(f"\n\n Velocity u value in the center of the first cell {u_cc[0,0,0]} ")
    
    # # print(f"velocities of the 1 vertex in the first cell: {u[0,0,0]:.6f}")
    # # print(f"velocities of the 2 vertex in the first cell: {u[1,0,0]:.6f}")
    # # print(f"velocities of the 3 vertex in the first cell: {u[1,0,1]:.6f}")
    # # print(f"velocities of the 4 vertex in the first cell: {u[1,1,1]:.6f}")
    # # print(f"velocities of the 5 vertex in the first cell: {u[1,1,0]:.6f}")
    # # print(f"velocities of the 6 vertex in the first cell: {u[0,1,1]:.6f}")
    # # print(f"velocities of the 7 vertex in the first cell: {u[0,1,0]:.6f}")
    # # print(f"velocities of the 8 vertex in the first cell: {u[0,0,1]:.6f}")
    
    # # print(f"Velocity average between the 8 vertex: {(u[0,0,0] + u[1,0,0] + u[1,0,1] +u[1,1,1] +u[1,1,0] +u[0,1,1] +u[0,1,0] +u[0,0,1])/8}")
    # print("======================================================================================================================================")
    
    dudxyz, dvdxyz , dwdxyz = compute_gradients_vel(u_cc,v_cc,w_cc,x_mid_vol,y_mid_vol,z_mid_vol)
    
    dudx, dudy, dudz = dudxyz
    dvdx, dvdy, dvdz = dvdxyz
    dwdx, dwdy, dwdz = dwdxyz
    
    # print("max value of gradient dudx: ",np.max(dudx))
    # print("max value of gradient dudy: ",np.max(dudy))
    # print("max value of gradient dudz: ",np.max(dudz))
    
    
    # print("max value of gradient dvdx: ",np.max(dvdx))
    # print("max value of gradient dvdy: ",np.max(dvdy))
    # print("max value of gradient dvdz: ",np.max(dvdz))
    
    
    # print("max value of gradient dwdx: ",np.max(dwdx))
    # print("max value of gradient dwdy: ",np.max(dwdy))
    # print("max value of gradient dwdz: ",np.max(dwdz))
    
    
    
    # print("\nmax value of velocity u_cc: ",np.max(u_cc))
    # print("max value of velocity v_cc: ",np.max(v_cc))
    # print("max value of velocity w_cc: ",np.max(w_cc))
    
    # print("\nmin value of dx :",np.min(dx_vol))
    # print("min value of dy :",np.min(dy_vol))
    # print("min value of dz :",np.min(dz_vol))
    
    
    
    # print(f"shapes of du/dx,dy,dz  :  {dudx.shape} {dudy.shape} {dudz.shape}" )
    
    Sxx = dudx[:-2,:-1,:-1]       ;        Syy = dvdy[:-2,:-1,:-1]      ;          Szz = dwdz[:-2,:-1,:-1]
    
    Sxy = 0.5 * (dudy[:-2,:-1,:-1] + dvdx[:-2,:-1,:-1])
    Sxz = 0.5 * (dudz[:-2,:-1,:-1] + dwdx[:-2,:-1,:-1])
    Syz = 0.5 * (dvdz[:-2,:-1,:-1] + dwdy[:-2,:-1,:-1])
    
    
    
    
    S_S = Sxx**2 + Syy**2 + Szz**2 + 2 * (Sxy**2 + Sxz**2 + Syz**2)
    vols = volumes[:-2,:-1,:-1]
    
    SS=integration_vol(S_S,vols)
    SR1 = SS*(2/1000)
    
    
    Re_uu = uu-(u*u) 
    Re_vv = vv-(v*v)
    Re_ww = ww-(w*w)
    Re_uv = uv-(u*v)
    Re_uw = uw-(u*w)
    Re_vw = vw-(vw*vw)
    
    Re_uu_cc,x_mid_vol, y_mid_vol, z_mid_vol=get_centered_cell_values(Re_uu, x, y, z)
    Re_vv_cc,x_mid_vol, y_mid_vol, z_mid_vol=get_centered_cell_values(Re_vv, x, y, z)
    Re_ww_cc,x_mid_vol, y_mid_vol, z_mid_vol=get_centered_cell_values(Re_ww, x, y, z)
    Re_uv_cc,x_mid_vol, y_mid_vol, z_mid_vol=get_centered_cell_values(Re_uv, x, y, z)
    Re_uw_cc,x_mid_vol, y_mid_vol, z_mid_vol=get_centered_cell_values(Re_uw, x, y, z)
    Re_vw_cc,x_mid_vol, y_mid_vol, z_mid_vol=get_centered_cell_values(Re_vw, x, y, z)
    
                                 #X direction                    #Y direction                              #Z direction      
    fst_term = Re_uu_cc[:-2,:-1,:-1] * dudx[:-2,:-1,:-1] + Re_uv_cc[:-2,:-1,:-1]*dvdx[:-2,:-1,:-1]+Re_uw_cc[:-2,:-1,:-1]*dwdx[:-2,:-1,:-1]
    snd_term = Re_uv_cc[:-2,:-1,:-1] * dudy[:-2,:-1,:-1] + Re_vv_cc[:-2,:-1,:-1]*dvdy[:-2,:-1,:-1]+Re_vw_cc[:-2,:-1,:-1]*dwdy[:-2,:-1,:-1]
    thrd_term= Re_uw_cc[:-2,:-1,:-1] * dudz[:-2,:-1,:-1] + Re_vw_cc[:-2,:-1,:-1]*dvdz[:-2,:-1,:-1]+Re_ww_cc[:-2,:-1,:-1]*dwdz[:-2,:-1,:-1]
    
    
    
    
    
    contraction=fst_term+snd_term+thrd_term
    
    re_stress = -1 * contraction
    total_Re_stress = integration_vol(contraction, vols)
    
   
    
    mu = 2/1000    
    tau_xx = 2 * mu * dudx
    tau_yy = 2 * mu * dvdy
    tau_zz = 2 * mu * dwdz
    tau_xy = mu * (dudy + dvdx)
    tau_xz = mu * (dudz + dwdx)
    tau_yz = mu * (dvdz + dwdy)
    
    #pause()
    
    E_shear= (-tau_xx*u[:-1,:-1,:-1]) + (-tau_xy*v[:-1,:-1,:-1]) + (-tau_xz*w[:-1,:-1,:-1])
    
    yz_areas,dy,dz=get_areas(z,y)
    Etau_cv=integration_yz(E_shear)
    Inlet_Outlet_Etau=Etau_cv[-1]#-Etau_cv[0]+Etau_cv[-1]
    
    xz_areas,dy,dz=get_areas(z,x)
    #Etau_cv_tb=integration_xz(E_shear)
    mult_areas = E_shear[:-1,:,:] * xz_areas[:-1,np.newaxis,:]
    sum_along_z=np.sum(mult_areas[:,:,:],axis=2)
    sum_along_x=np.sum(sum_along_z[:,:],axis=0)
    Etau_cv_tb=sum_along_x
    Top_bottom_Etau=-Etau_cv_tb[-1]+Etau_cv_tb[0]
    
    
    #pause()
    
    arf_u = u_fc_reshape-1.00037          ;      arf_v = v_fc_reshape-0.00021       ;       arf_w = w_fc_reshape-0    ;      p_inf=0
    
    
    
    
    Ek = 0.5*(arf_u**2)*u_fc_reshape #Ea
    Ev = 0.5*(arf_v**2+arf_w**2)*u_fc_reshape
    Ep = (  (p_fc_reshape)    +p_inf)*arf_u
    
    #print(Ep)
    
    yz_areas,dy,dz=get_areas(z,y)
    Ek_cv=integration_yz(Ek)
    Ev_cv=integration_yz(Ev)
    
    yz_areas_na=yz_areas[np.newaxis,:,:]
    mult_areas = Ep * yz_areas_na
    
    sum_along_y=np.sum(mult_areas[:,:,:],axis=-1)
    
    sum_along_z=np.sum(sum_along_y[:,:],axis=-1)
    print(sum_along_z)
    
    Ep_cv=sum_along_z ;
    
    
    
    Ep_cv=integration_yz(Ep)
    
    #print(Ep_cv)
    
    xz_areas,dy,dz=get_areas(z,x)
    Ek_cv_tb=integration_xz(Ek)
    Ev_cv_tb=integration_xz(Ev)
    Ep_cv_tb=integration_xz(Ep)
    
    Inlet_Outlet_Ek= Ek_cv[-1]#-Ek_cv[0]+Ek_cv[-1]
    Inlet_Outlet_Ev= Ev_cv[-1]#-Ev_cv[0]+Ev_cv[-1]
    Inlet_Outlet_Ep= Ep_cv[-1]#-Ep_cv[0]+Ep_cv[-1]
    Top_bottom_Ek= 0*(-Ek_cv_tb[-1]+Ek_cv_tb[0])
    Top_bottom_Ev= 0*(-Ev_cv_tb[-1]+Ev_cv_tb[0])
    Top_bottom_Ep= 0*(-Ep_cv_tb[-1]+Ep_cv_tb[0])
    
    Inlet_Outlet_Ek_perc = (-Ek_cv[0]+Ek_cv[-1])/-Ek_cv[0]
    Inlet_Outlet_Ev_perc = (-Ev_cv[0]+Ev_cv[-1])/-Ev_cv[0]
    Inlet_Outlet_Ep_perc = (-Ep_cv[0]+Ep_cv[-1])/-Ep_cv[0]
    Top_bottom_Ek_perc = (-Ek_cv_tb[-1]+Ek_cv_tb[0])/Ek_cv_tb[-1]
    Top_bottom_Ev_perc = (-Ev_cv_tb[-1]+Ev_cv_tb[0])/Ev_cv_tb[-1]
    Top_bottom_Ep_perc = (-Ep_cv_tb[-1]+Ep_cv_tb[0])/Ep_cv_tb[-1]
    
    
    
    
    total_anergy = total_Re_stress+SR1
    
    Ea_total = Inlet_Outlet_Ek+Top_bottom_Ek
    Ev_total = Inlet_Outlet_Ev+Top_bottom_Ev
    Ep_total = Inlet_Outlet_Ep+Top_bottom_Ep
    Etau_total=Inlet_Outlet_Etau+Top_bottom_Etau
    total_mech_exergy = Ea_total+Ev_total+Ep_total+Etau_total
    
    
    Force_normalisation = 0.5*0.1666667
    
    print("\n\n\n====================================================================================================================================")
    print("Mechanical Exergy Values")
    print("====================================================================================================================================")
    # print(f"\nKinetic Energy lost or gained in the control volume:    {(Inlet_Outlet_Ek+Top_bottom_Ek):.4f} , { ( (Inlet_Outlet_Ek_perc+Top_bottom_Ek_perc)*100  ):.4f} %")
    # print(f"\nTransverse Energy lost or gained in the control volume: {(Inlet_Outlet_Ev+Top_bottom_Ev):.4f} , { ( (Inlet_Outlet_Ev_perc+Top_bottom_Ev_perc)*100  ):.4f} %")
    # print(f"\nPressure Energy lost or gained in the control volume:   {(Inlet_Outlet_Ep+Top_bottom_Ep):.4f} , { ( (Inlet_Outlet_Ep_perc+Top_bottom_Ep_perc)*100  ):.4f} %")
    print("\nLoss of energy due to mean kinetic energy: ",total_Re_stress)
    print("\nTotal viscous dissipation rate: ",total_anergy)
    print("\nTotal mechanical exergy : ",total_mech_exergy)
    print("\n\n\nCd due to mechanical exergy: ",(total_mech_exergy)/(Force_normalisation))
    
    print("\nEa :",Ea_total)
    print("Ev :",Ev_total)
    print("Ep :",Ep_total)
    print("Ek :",Ea_total+Ev_total)
    
    
    print("\nCd due to anergy: ",(total_anergy)/(Force_normalisation))
    print("Cd due to exergy: ",(total_mech_exergy)/(Force_normalisation))
    print("Total Cd : ", (total_anergy+total_mech_exergy)/(Force_normalisation))
    print("====================================================================================================================================\n\n\n\n")
    
#    input("Press Enter to continue...")


    titles=["Ea","Ev","Ep","Etau","Ek","Anergy","Exergy","Cd total"]
    
    Ea_plot.append(Ea_total/Force_normalisation)
    Ev_plot.append(Ev_total/Force_normalisation)
    Ep_plot.append(Ep_total/Force_normalisation)
    Etau_plot.append(Etau_total/Force_normalisation)
    Ek_plot.append((Ea_total + Ev_total)/Force_normalisation)
    Anergy_plot.append(total_anergy/Force_normalisation)
    Exergy_plot.append(total_mech_exergy/Force_normalisation)
    Cd_tot_plot.append((total_anergy + total_mech_exergy)/Force_normalisation)
    x_coords_plot.append(cv_x_max)
    print("integrating variables at : ",cv_x_max)
    
    
 

with open('Anergy_Exergy.dat', mode='w') as file:
    # Write the header for all columns
    
    header = ['x-coords'] + titles
    file.write('\t\t '.join(header) + '\n')


    # Write the data rows
    for i in range(len(x_coords_plot)):
        row = [x_coords_plot[i], Ea_plot[i], Ev_plot[i], Ep_plot[i], Etau_plot[i], Ek_plot[i], Anergy_plot[i], Exergy_plot[i], Cd_tot_plot[i]]
        file.write('\t\t'.join(map(str, row)) + '\n')
        

  

    
plt.plot( x_coords_plot  , Ea_plot, label="(Ea))")
plt.plot( x_coords_plot  , Ev_plot, label="(Ev))")
plt.plot( x_coords_plot  , Ep_plot, label="(Ep))")
plt.plot( x_coords_plot  , Etau_plot, label="(Etau))")
plt.plot( x_coords_plot  , Ek_plot, label="(Ek))")
plt.plot( x_coords_plot  , Exergy_plot, label="Exergy (E))")
plt.plot( x_coords_plot  , Anergy_plot, label="Anergy (phi))")
plt.plot( x_coords_plot  , Cd_tot_plot, label="Total Cd")

plt.axvline(x=cv_x_min+(2*prm.dx), color='gray', linestyle=':', linewidth=1, alpha=0.9)
plt.axvline(x=cv_x_min, color='red', linestyle=':', linewidth=1, alpha=0.9)
plt.axvline(x=cv_x_max_0, color='red', linestyle=':', linewidth=1, alpha=0.9, label="Airfoil trailing edge")
plt.axvline(x=cv_x_max, color='gray', linestyle=':', linewidth=1, alpha=0.9)
plt.axhline(y=0.1438,color="blue",linestyle='--',linewidth=1,alpha=0.9, label = "Cd in RRF")
plt.axhline(y=0,color="black",linestyle='--',linewidth=1,alpha=0.9, label = "Cd in ARF")

# #plt.plot( x, (SR_along_x-Re_along_x)/(0.5*0.1666667) ,label="phi"   )
plt.xlabel("x-coordinates")
plt.ylabel("Cd")
plt.title("Cd vs x-coordinates")
plt.legend(framealpha=0.0)
plt.minorticks_on()
plt.grid(which='both', color='gray', linestyle='-', linewidth=0.5, alpha=0.2)
plt.tight_layout()
plt.legend()
plt.savefig("Energy_plot.png")
plt.show()
pause()