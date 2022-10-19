# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 22:33:18 2022.

@author: RUBIN

init file module for MGD.py
"""
from numpy import random

# %% Parameters
# =============================================================================
# Set Initial parameters (please modify here if necessary)
# =============================================================================
domain_shape_profile = 'Square'
characteristic_length = 1

x_initial_coordinate = 0.0 #Axial direction
y_initial_coordinate = 0.0 #Vertical direction
z_initial_coordinate = 0.0 #Third direction

time_initial_value = 0.0
timestep_size = 0.01
max_timesteps = 100

initial_velocity_distribution = 'Normal' #Normal/Random
max_initial_velocity = 100

wall_reflection_scheme = 'Specular' #Specular/Diffuse

# %% Shape Function
# =============================================================================
# Define Shape Function 
# =============================================================================
def domain_shape(shape):
    """
    
    Define shape of domain. Returns integer to represent chosen shape.
    
          Y

          |
          |
          |________  X
          /
         /
        /
      Z 
      
    Parameters
    ----------
    shape : string
        Chosse between 'Circle', 'Square' and 'Hex'.

    Returns
    -------
    shape_func : int
        0,1,2 are circular, square and hex-cross sectioned pipe respectively.

    """
    if shape == 'Circle': #Pipe with circular cross-section
        shape_func = 0
    elif shape == 'Square': #Pipe with square cross-section
        shape_func = 1
    elif shape == 'Hex':    #Pipe with hexagonal cross-section
        shape_func = 2
    return shape_func

shape_function = domain_shape(domain_shape_profile)

# %% Useful Functions
# =============================================================================
# Function Definitiions
# =============================================================================
def diff_vel_def(vel_dist, max_vel, init_flag):
    """
    Set velocities for diffuse velocity distribution at reflection.

    Parameters
    ----------
    vel_dist : string
        choice of velocity distribution - Normal(Gaussian)/Random.
    max_vel : float
        maximum velocity used for calculation.
    init_flag : int
        flag = 1 when describing initial velocity dist., flag=0 otherwise.

    Returns
    -------
    xvel : float
        output velocity in x-direction.
    yvel : float
        output velocity in y-direction.
    zvel : float
        output velocity in z-direction.

    """
    if (vel_dist == 'Normal'):
        xvel = random.normal(loc = 0.0, scale = max_vel/3)
        yvel = random.normal(loc = 0.0, scale = max_vel/3)
        zvel = random.normal(loc = 0.0, scale = max_vel/3)
    elif (vel_dist == 'Random'):
        xvel = random.uniform(-max_vel, max_vel)      
        yvel = random.uniform(-max_vel, max_vel)      
        zvel = random.uniform(-max_vel, max_vel)      
    if (init_flag == 1):
        xvel = abs(xvel)
    return xvel, yvel, zvel

def spec_vel_def(u, v, w, normal):
    """
    Define output velocities after specular reflection.

    Parameters
    ----------
    u : float
        velocity in x-direction.
    v : float
        velocity in y-direction.
    w : float
        velocity in z-direction.
    normal : list/array [x, y, z]
        normal vector direction cosines.

    Returns
    -------
    u_out : float
        reflected x-direction velocity.
    v_out : float
        reflected y-direction velocity.
    w_out : float
        reflected z-direction velocity.

    """
    from MGD_init import shape_function
    u_out = u
    if (shape_function == 0): #Circle case
        pass
    elif (shape_function == 1): #Square case
        if (normal[2] == 0): #normal in vertical direction ie top/bot wall
            v_out = -v
            w_out = w
        elif (normal[1] == 0): #normal in horiz direction i.e side wall
            v_out = v
            w_out = -w
    elif (shape_function == 2): #Hex case
        pass
    #v_out = 0 #check
    #w_out = 0 #check
    return u_out, v_out, w_out

def sign(input_var):
    """
    User-defined Signum function. Outputs array if array inputted.

    Parameters
    ----------
    input_var : float/float array
        innput to Signum function.

    Returns
    -------
    output_var : int/int array
        Output of Signum function.

    """
    output_var = input_var
    if hasattr(input_var, "__len__"):#if array/list
        for i in range(len(input_var)):
            if (input_var[i] == 0.0):
                output_var[i] = 0
            elif (input_var[i] > 0.0):
                output_var[i] = 1
            elif (input_var[i] < 0.0):
                output_var[i] = -1
    
    else:#if NOT array/list
        if (input_var == 0.0):
            output_var = 0
        elif (input_var > 0.0):
            output_var = 1
        elif (input_var < 0.0):
            output_var = -1
    return output_var
    
def zero_correct(input_var):
    """
    Replace zero with small value to prevent division by zero.

    Parameters
    ----------
    input_var : int/float/float-array
        input variable/array.

    Returns
    -------
    output_var : same type as input_variable (int becomes float)
        corrected output for input=zero.

    """
    error_tolerance = 1.0e-16
    output_var = input_var
    if hasattr(input_var, "__len__"): #if array
        for i in range(len(input_var)):
            if (input_var[i] == 0.0):
                output_var[i] = error_tolerance
    else: #if NOT array
        if (input_var == 0.0):
                output_var = error_tolerance
    return output_var

def reflection_select(scheme, velx, vely, velz, normal):
    """
    Call appropriate reflection scheme. Conceived to create cleaner code.

    Parameters
    ----------
    scheme : string
        'Specular' to call spec_vel_def function
        'Diffuse' to call diff_vel_def function
    velx : float
        velocity in axial direction
    vely : float
        velocity in vertical direction
    velz : float
        velocity in third direction
    normal: float array
        wall normal vector expressed as its direction cosines [nx, ny, nz]

    Returns
    -------
    u_out : float
        post-reflection velocity in axial direction (x)
    v_out : float
        post-reflection velocity in vertical direction (y)
    w_out : float
        post-reflection velocity in third direction (z)
    """
    if (scheme == 'Specular'):
        u_out, v_out, w_out = spec_vel_def(velx, vely, velz, normal)
    elif (scheme == 'Diffuse'):
        u_out, v_out, w_out = diff_vel_def(initial_velocity_distribution, \
                                           max_initial_velocity, 0)
    return u_out, v_out, w_out

# %% Initialization    
# =============================================================================
# Define inlet velocities
# =============================================================================
u_inlet, v_inlet, w_inlet = diff_vel_def('Normal', max_initial_velocity, 1)