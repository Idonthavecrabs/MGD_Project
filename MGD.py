# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 14:59:41 2022.

@author: Rubin

@version: 1.0

Molecular gas dynamics (MGD):
Code to track a particle's trajectory through domain of specified shape
"""

import numpy as np
import random
from datetime import datetime
from MGD_init import *



# %% Built-in Functions
# =============================================================================
# Function Definitions
# =============================================================================

def collision(x_old,y_old,z_old,x_vel,y_vel,z_vel,ts,coll,repeat):
    """
    Check if collision with wall exists. If yes, then also returns point of collision with wall with time-correction for calculating trajectory of post-collision.
    
    Parameters
    ----------
    x_old : old-x position value (float)
        .
    y_old : old-y position value (float)
        .
    z_old : old-z position value (float)
        .
    x_vel : x-direction velocity value (float)
        Value is pre-correction while checking for wall collisions.
    y_vel : y-direction velocity value (float)
        Value is pre-correction while checking for wall collisions.
    z_vel : z-direction velocity value (float)
        Value is pre-correction while checking for wall collisions.
    coll  : list/array
        int-Array that keeps track of collisions throughout 0=no collision
        1=collision, Note that each entry in array represents an instance.
        Number of instances >= Number of timesteps. A single timestep may
        have multiple collisions i.e multiple instances in that timestep
    ts    : timesize (float)
        Total time for calculating final position = timestep if no collision
        ts = timestep - (time taken for previous collisions in current 
                         timestep)  ==> if called recursively
    ##shape_func : string
        ##Describe domain cross-section shape.
    ##d     : characteristic length (float)
        ##Characteristic length of domain: dia for circle, side length for
        ##square, side length for hexagon.
    repeat : int
        Flag for recursive call. 0 = First function call
                                 1 = Recursive call within function

    Internal Variables
    ------------------
    t_coll : (float) time taken for collision with wall during corrent timestep
        if applicable
    t_corr : (float) corrected time = timestep - t_coll
    normal : float array = [n_x,n_y,n_z]
        array containing direction cosines of normal at point of collision

    Returns
    -------
    None.

    """
    from MGD_init import wall_reflection_scheme, max_initial_velocity, \
                            initial_velocity_distribution, spec_vel_def, \
                            diff_vel_def, zero_correct, sign, \
                            reflection_select, shape_function, \
                            characteristic_length
    #from MGD_init import *
    
    d = characteristic_length
    if (repeat == 0):
        wall = 0
        dist = -1.0
        
    x_corr, y_corr, z_corr = x_old + x_vel*ts, y_old + y_vel*ts, \
                             z_old + z_vel*ts #Corrected Positions
    t_corr = 0 #Time correction
    x_vel = zero_correct(x_vel)
    y_vel = zero_correct(y_vel)
    z_vel = zero_correct(z_vel)
    #Temporary initialization
    u_new = x_vel
    v_new = y_vel
    w_new = z_vel
    #Cases for profile shape:
    while (ts > 0.0):
        if shape_function == 0: #Circle Pipe
            measure = y_corr**2 + z_corr**2 - d**2
            if (measure >= 0): #condition for isinsidewall fails at >=0
                wall = 1
                coll.append(1) #register collision in coll array
                t_coll = 1/(y_vel**2 + z_vel**2)*(-(y_vel*y_old + z_vel*z_old)\
                         + ((y_vel*y_old + z_vel*z_old)**2 + ((d/2)**2 \
                         - y_old**2 - z_old**2)*(y_vel**2 + z_vel**2))**0.5)
                t_corr = ts - t_coll
                x_coll, y_coll, z_coll = x_old + x_vel*t_coll, y_old + \
                                        y_vel*t_coll, z_old + z_vel*t_coll
               # normal is given by y/y1 = z/z1 where z1,y1 are point on circle 
               #from where normal is drawn
                normal = [0, y_coll*2/d, z_coll*2/d]
                if (wall_reflection_scheme == 'Specular'):
                    u_new, v_new, w_new = reflection_select('Specular',x_vel, \
                                                        y_vel, z_vel, normal)
                elif (wall_reflection_scheme == 'Diffuse'):
                    u_new, v_new, w_new = reflection_select('Diffuse',x_vel, \
                                                        y_vel, z_vel, normal)
                x_corr, y_corr, z_corr = x_coll + u_new*t_corr, y_coll + \
                                        v_new*t_corr, z_coll + w_new*t_corr
                collision(x_coll,y_coll,z_coll,u_new,v_new,w_new,t_corr,coll,1)                        
                
        elif shape_function == 1: #Square Pipe
            measure = max(abs(y_corr) - d/2, abs(z_corr) - d/2)
            if (measure >= 0):
                wall = 1
                coll.append(1)
                t_coll_y = abs(d/2 - sign(y_vel)*y_old)/y_vel
                t_coll_z = abs(d/2 - sign(z_vel)*z_old)/z_vel
                if (abs(t_coll_y) < abs(t_coll_z)):
                    t_coll = abs(t_coll_y)
                    normal = [0, sign(t_coll_y), 0]
                else:
                    t_coll = abs(t_coll_z)
                    normal = [0, 0, sign(t_coll_z)]
                t_corr = ts - t_coll
                x_coll, y_coll, z_coll = x_old + x_vel*t_coll, y_old + \
                                        y_vel*t_coll, z_old + z_vel*t_coll
                if (wall_reflection_scheme == 'Specular'):
                    u_new, v_new, w_new = reflection_select('Specular',x_vel, \
                                                        y_vel, z_vel, normal)
                elif (wall_reflection_scheme == 'Diffuse'):
                    u_new, v_new, w_new = reflection_select('Diffuse',x_vel, \
                                                        y_vel, z_vel, normal)
                x_corr, y_corr, z_corr = x_coll + u_new*t_corr, y_coll + \
                                        v_new*t_corr, z_coll + w_new*t_corr
                collision(x_coll,y_coll,z_coll,u_new,v_new,w_new,t_corr,coll,1)                        
    
        elif shape_function == 2: #Hex Pipe
            pass
    
    if (wall == 0):
        coll.append(0)

    return x_corr, y_corr, z_corr, x_coll, y_coll, z_coll, u_new, v_new, \
            w_new, t_corr, coll

def second_func():
    pass



# %% Initialization 

now = datetime.now()
print("Start time:", end='')
print(now.strftime("%d-%m-%Y  %H:%M:%S:%f"))

# =============================================================================
# Initialize - Take values from MGD_init file
# =============================================================================
x_init, y_init, z_init, time_init = x_initial_coordinate, y_initial_coordinate\
                                    , z_initial_coordinate, time_initial_value
timestep = timestep_size
x, y, z, coll= [], [], [], []
time = 0
time_final = max_timesteps * timestep
x.append(x_init)
y.append(y_init)
z.append(z_init)
u_init = zero_correct(u_inlet)
v_init = zero_correct(v_inlet)
w_init = zero_correct(w_inlet)

chosen_shape = domain_shape_profile
d = characteristic_length 

# %% Timestepping
# =============================================================================
# Step through timesteps
# =============================================================================
while (time <= time_final):
    time += timestep
    

