"""
	plot_thrust_profile.py

	Purpose: Implements the functions to produce thrust profiles of tranfers.

	@author Thomas Caleb

	@version 1.0 16/01/2024
    
"""

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from classes import Dataset

"""
    Plots a thrust profile with a stairs shape.

"""
def plot_stairs_profile(dt, u, ax, label, color,
                        linestyle="solid"):
    nb_point_stair = 50
    nb_point = len(dt)
    t = []
    u_stairs = []
    current_time = 0
    for i in range(nb_point):
        dt_i = dt[i]/(float(nb_point_stair))
        u_i = u[i]
        for j in range(nb_point_stair):
            t.append(current_time)
            u_stairs.append(u_i)
            current_time = current_time + dt_i
            
    ax.plot(t, u_stairs, label=label,
            color=color, linestyle=linestyle)

"""
    Plots a thrust pofile for a given transfer dataset.

"""
def plot_thrust_profile(dataset):
    # Settings
    
    dpi = 200
    
    # Thrust norm
    color_thrust_norm = "black"
    label_thrust_norm  = "Thrust"
    
    # Thrust components
    plot_components = False
    color_thrust_0 = "green"
    color_thrust_1 = "cyan"
    color_thrust_2 = "grey"
     
    # Thrust norm
    color_thrust_max = "red"
    linestyle_thrust_max = "dotted"
    label_thrust_max  = "Max thrust"
    
    # Normalisation
    denormalise = True
    
    # Legend
    show_legend = False
    legend_loc = "lower right"

    # Output
    show_grid = True
    save_figure = True
    saving_format = "pdf"
    show_plot = True
    
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Control"):
            data_control = dataset.list_datasets[i].copy()
    ux = data_control[0,:]
    uy = data_control[1,:]
    uz = data_control[2,:]
    dt = data_state[7,:-1]
    u = np.sqrt(ux*ux + uy*uy + uz*uz)
    max_trust = dataset.spacecraft_parameters.thrust
    x_label = "Time [TU]"
    y_label = "Thrust norm [THRUSTU]"
    
    # Normalisation
    if denormalise:
        THRUSTU = dataset.spacecraft_parameters.constants.THRUSTU
        TU = dataset.spacecraft_parameters.constants.TU
        ux *= THRUSTU
        uy *= THRUSTU
        uz *= THRUSTU
        u *= THRUSTU
        max_trust *= THRUSTU
        dt *= (TU/86400)
        
        # Labels
        x_label = x_label.replace("TU", "days")
        y_label = y_label.replace("THRUSTU", "N")


    # Create plot
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot()

    # Set labels
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    
    # Axies limites
    plt.ylim((0 - 0.05*max_trust, 1.05*max_trust))
    
    # Plot Thrust 
    plot_stairs_profile(dt, u, ax, label_thrust_norm, color_thrust_norm)
    
    # Separate components
    if plot_components:
        plot_stairs_profile(dt, ux, ax, "ux", color_thrust_0)
        plot_stairs_profile(dt, uy, ax, "uy", color_thrust_1)
        plot_stairs_profile(dt, uz, ax, "uz", color_thrust_2)
    
    # Plot max thrust
    plot_stairs_profile(dt, u*0 + max_trust, ax,
                        label_thrust_max, color_thrust_max,
                        linestyle_thrust_max)
    
    if show_grid:
        plt.grid()
    
    if show_legend: 
        plt.legend(loc=legend_loc)
    
    if save_figure:
        # Get adress
        file_name = dataset.file_name.replace(
            "datasets", "plots")
                
        # Add signature
        signature = ("_thrust" + "." + saving_format)
        file_name = file_name.replace(
            ".dat", signature)
        
        # Save
        plt.savefig(file_name)    
       
    if show_plot:   
        plt.show()
    else:
        plt.close(fig)
    