#! python3 
# *****************************************************************************
#
# File:                 streamlinesForMOPS.py
# Project:              python3
# Author(s):            Eric J. Bringley (eb656)
#
# Contact:
#   Mr. Eric J. Bringley
#   Department of Chemical Engineering
#   University of Cambridge
#   West Cambridge Site
#   Philippa Fawcett Drive
#   Cambridge
#   CB3 0AS
#   United Kingdom
#
#   Email:   eb656@cam.ac.uk
#   Website: como.cheng.cam.ac.uk
#
# Purpose:
#   This script takes streamlines from OpenFOAM and formats them as inputs for
#   the TiO2 (TTIP) particle models in MOPS.
#
# Questions for Casper/Noel:
#   1. Velocities are reported in magnitudes, therefore, I assume this is the 
#       magnitude of the velocity projected into the streamline. Is this correct?
#   2. Time is reported in two locations: Mops.inx and gasphase.csv -- and they 
#       have the same values... Does this include U and U_thermophoretic?
#
# *****************************************************************************

import sys, os, re
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
# import xml.etree.ElementTree as ET

# Functions to assist calculations
def project(vector1, vector2):
    "Project vector1 onto vector 2"
    # check lengths of vectors:
    vector2_length = np.sqrt(dot(vector2,vector2))
    vector2_unit = vector2/vector2_length
    projection = dot(vector1,vector2_unit)
    return projection

def dot(v1,v2):
    if (len(v1)!= 2):   
        raise ValueError
    if (len(v2)!= 2):   
        raise ValueError
    return v1[0]*v2[0] + v1[1]*v2[1]


# Inputs: 
foam_case_dir = os.path.abspath("/rds/user/eb656/hpc-work/Stagnation2/20200619_mym_phi1.0/")
foam_time = 0.2
foam_prostProcess_dir = "postProcessing/sets/streamLines/"
foam_streamline_basename = "track0"
foam_streamline_dir = os.path.join(foam_case_dir,foam_prostProcess_dir+str(foam_time))
foam_streamline_file_scalar = os.path.join(foam_streamline_dir, foam_streamline_basename+"_p_T.csv" )
foam_streamline_file_vector = os.path.join(foam_streamline_dir, foam_streamline_basename+"_U_V_T.csv" )


# Read Streamine Velocity Data
velocity_data = pd.read_csv(foam_streamline_file_vector)
velocity_data = velocity_data.drop(columns=['x'])
scalar_data = pd.read_csv(foam_streamline_file_scalar)
scalar_data = scalar_data.drop(columns=['x'])

### separate streamlines:
# Assumptions:
#   1. Streamlines start at their max height and travel downward
#   2. Streamlines do not drift to above the starting height

## Check raw streamlines as XY plot: 
checkRawStreamlines = False
if(checkRawStreamlines):
    velocity_data.plot(x='y', y='z', kind='line')
    plt.show()

# find height of starting streamline 
start_z = velocity_data['z'].max()
print("Starting height {}".format(start_z))
streamline_idx = np.empty_like(velocity_data['z'])

# Calculate streamline indexes:
first_steamline_index=0
idx=first_steamline_index-1

# streamline info: 
t_res = 0
u_mag = 0

streamline_y = velocity_data['y']
streamline_z = velocity_data['z']
streamline_uy = velocity_data['U_1'] # velocity y
streamline_uz = velocity_data['U_2'] # velocity z
streamline_vy = velocity_data['V_T_1'] # thermophoretic velocity y
streamline_vz = velocity_data['V_T_2'] # thermophoretic velocity z

for i in range(len(streamline_z)):
    # Strreamline idx:
    height = streamline_z[i]
    if (height == start_z):
        print("New streamline")
        idx = idx+1
    streamline_idx[i]=idx
    # print("Height is {} in steamline {}".format(height,idx))

#
velocity_data['streamline'] = streamline_idx
velocity_data['streamline'] = velocity_data['streamline'].astype('int') # convert to integer!
# velocity_data['streamline'] = streamline_idx
scalar_data['streamline']   = streamline_idx
scalar_data['streamline'] = scalar_data['streamline'].astype('int')

# Calculate residence time and velocity magnitudes:
umag_full = np.zeros(0)
vmag_full = np.zeros(0)
tres_full = np.zeros(0)

# direction is calculated as change from before and after points.
# distance is change from previous point to current point

for j in range(0,velocity_data['streamline'].astype('int').max()+1):
    print("Post Processing Streamline {}".format(j))
    dfStreamline = velocity_data[velocity_data['streamline'] == j]
    streamline_y = dfStreamline['y'].to_numpy()
    streamline_z = dfStreamline['z'].to_numpy()
    streamline_uy = dfStreamline['U_1'].to_numpy()# velocity y
    streamline_uz = dfStreamline['U_2'].to_numpy() # velocity z
    streamline_vy = dfStreamline['V_T_1'].to_numpy() # thermophoretic velocity y
    streamline_vz = dfStreamline['V_T_2'].to_numpy() # thermophoretic velocity z
    streamline_umag = np.empty_like(dfStreamline['z'])
    streamline_vmag = np.empty_like(dfStreamline['z'])
    streamline_tres = np.empty_like(dfStreamline['z'])
    # intial point:
    position_i = [streamline_y[0], streamline_z[0]]
    u = [streamline_uy[0], streamline_uz[0]]
    v = [streamline_vy[0], streamline_vz[0]]
    position_next = [streamline_y[1], streamline_z[1]]
    on_to_direction = [position_next[0]-position_i[0], position_next[1]-position_i[1]]
    streamline_umag[0] = project(u,on_to_direction)
    streamline_vmag[0] = project(v,on_to_direction)
    streamline_tres[0] = 0
    # internal points:
    for i in range(1,len(streamline_z)-1):
        position_i = [streamline_y[i], streamline_z[i]]
        u = [streamline_uy[i], streamline_uz[i]]
        v = [streamline_vy[i], streamline_vz[i]]
        position_prev = [streamline_y[i-1], streamline_z[i-1]]
        position_next = [streamline_y[i+1], streamline_z[i+1]]
        on_to_direction = [position_next[0]-position_prev[0], position_next[1]-position_prev[1]]
        u_v = [u[0]+v[0], u[1]+v[1]]
        u_total = project(u_v,u_v)
        dx = [position_i[0]-position_prev[0], position_i[1]-position_prev[1]]
        distance = project(dx,dx)
        tres_i = distance/u_total
        streamline_umag[i] = project(u,on_to_direction)
        streamline_vmag[i] = project(v,on_to_direction)
        streamline_tres[i] = tres_i+streamline_tres[i-1]
    # final point:
    position_i = [streamline_y[-1], streamline_z[-1]]
    position_prev = [streamline_y[-2], streamline_z[-2]]
    u = [streamline_uy[-1], streamline_uz[-1]]
    v = [streamline_vy[-1], streamline_vz[-1]]
    on_to_direction = [position_i[0]-position_prev[0], position_i[1]-position_prev[1]]
    u_v = [u[0]+v[0], u[1]+v[1]]
    u_total = project(u_v,u_v)
    dx = [position_i[0]-position_prev[0], position_i[1]-position_prev[1]]
    distance = project(dx,dx)
    tres_i = distance/u_total
    streamline_umag[-1] = project(u,on_to_direction)
    streamline_vmag[-1] = project(v,on_to_direction)
    streamline_tres[-1] = streamline_tres[-2]+ tres_i
    # print(streamline_umag)
    # append data to full lists:
    umag_full = np.concatenate([umag_full, streamline_umag])
    vmag_full = np.concatenate([vmag_full, + streamline_vmag])
    tres_full = np.concatenate([tres_full, + streamline_tres])

print(len(umag_full))
print(len(streamline_idx))
velocity_data['ConvectiveVelocity[m/s]'] = umag_full
velocity_data['ThermophoreticVelocity[m/s]'] = vmag_full
velocity_data['Time'] = tres_full

print(velocity_data)
print(velocity_data.dtypes)

# Create master dataframe: