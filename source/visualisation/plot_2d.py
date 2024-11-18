"""
	plot_2d.py

	Purpose: Implements the functions to produce 2D plots of tranfers.

	@author Thomas Caleb

	@version 1.0 16/01/2024
    
"""

import numpy as np
import matplotlib.pyplot as plt
#plt.switch_backend('agg')
import scipy.interpolate as interpolate

from misc import get_Lagrange_point
from classes import Dataset


"""
    Plots the specific points for a given dynamical system (two-body and 
    circular-restricted three-body problem).

"""
def plot_system_points(dataset, axis_0, axis_1, ax,
                       list_colors, list_markers,
                       list_plots, denormalise):    
    # Discriminate dynamical system
    if dataset.dynamical_system.startswith("TBP"):
        # Make points
        coord_center = np.array([0, 0, 0])
        
        # Denormalise
        if denormalise:
            LU = dataset.spacecraft_parameters.constants.LU
            coord_center *= LU
        
        # Plot
        if list_plots[0]: # Central body
            ax.scatter(coord_center[axis_0], coord_center[axis_1],
                       color=list_colors[0], marker=list_markers[0])
            
            # Print name
            central_body_name = dataset.dynamical_system.split(" ")[1]
            ax.text(coord_center[axis_0], coord_center[axis_1],
                    " " + central_body_name)

    elif dataset.dynamical_system.startswith("CR3BP"):
        # Make points
        coord_P1 = np.array([-dataset.spacecraft_parameters.constants.MU, 0, 0])
        coord_P2 = np.array([1-dataset.spacecraft_parameters.constants.MU, 0, 0])
        coord_L1 = get_Lagrange_point(dataset, 1)
        coord_L2 = get_Lagrange_point(dataset, 2)
        
        # Denormalise
        if denormalise:
            LU = dataset.spacecraft_parameters.constants.LU
            coord_P1 *= LU
            coord_P2 *= LU
            coord_L1 *= LU
            coord_L2 *= LU
            
        # Get names
        primary_names = dataset.dynamical_system.split(" ")[1].split("-")

        # Plots
        if list_plots[0]: # First primary
            ax.scatter(coord_P1[axis_0], coord_P1[axis_1],
                       color=list_colors[0], marker=list_markers[0])
            ax.text(coord_P1[axis_0], coord_P1[axis_1],
                    " " + primary_names[0])
            
        
        if list_plots[1]: # Second primary
            ax.scatter(coord_P2[axis_0], coord_P2[axis_1],
                       color=list_colors[1], marker=list_markers[1])
            ax.text(coord_P2[axis_0], coord_P2[axis_1],
                    " " + primary_names[1])
            
        if list_plots[2]: # Lagrange point 1
            ax.scatter(coord_L1[axis_0], coord_L1[axis_1],
                       color=list_colors[2], marker=list_markers[2])
            ax.text(coord_L1[axis_0]+0.02, coord_L1[axis_1]-0.075,
                    " $L_1$")
        
        if list_plots[3]: # Lagrange point 2
            ax.scatter(coord_L2[axis_0], coord_L2[axis_1],
                       color=list_colors[3], marker=list_markers[3])  
            ax.text(coord_L2[axis_0]-0.1, coord_L2[axis_1]-0.075,
                    " $L_2$")
        
"""
    Plots the departure and arrival points of a transfer.

"""    
def plot_departure_arrival(dataset, axis_0, axis_1, ax,
                           list_colors, list_markers,
                           denormalise):    
    # Retrieve data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i].copy()
            
    # Denormalise
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        data_state *= LU
            
    # Plot departure and arrival
    ax.scatter(data_state[axis_0, 0], data_state[axis_1, 0],
               color=list_colors[0], marker=list_markers[0])
    ax.text(data_state[axis_0, 0], data_state[axis_1, 0],
            " $x_0$")
    ax.scatter(data_state[axis_0, -1], data_state[axis_1, -1],
               color=list_colors[1], marker=list_markers[1])
    ax.text(data_state[axis_0, -1], data_state[axis_1, -1],
            " $x_t$")
    
"""
    Plots the thrust vectors along the trajectory.

"""
def plot_thrust_vector(dataset, axis_0, axis_1, ax,
                       thrust_color,
                       denormalise):
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Control"):
            data_control = dataset.list_datasets[i].copy()
    coord_0 = data_state[axis_0,:]
    coord_1 = data_state[axis_1,:]
    ucoord_0 = data_control[axis_0,:]
    ucoord_1 = data_control[axis_1,:]
    
    # Denormalise
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        THRUSTU = dataset.spacecraft_parameters.constants.THRUSTU
        coord_0 *= LU
        coord_1 *= LU
        ucoord_0 *= THRUSTU
        ucoord_1 *= THRUSTU
    
    # Plot arrows
    ax.quiver(coord_0[:-1], coord_1[:-1],
              ucoord_0, ucoord_1,
              color=thrust_color, label='Thrust')
    
"""
    Plots the departure and arrival orbits of a transfer.

"""
def plot_reference_orbits(dataset, axis_0, axis_1, ax,
                          list_colors, list_linestyles,
                          denormalise):    
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Departure orbit"):
            data_departure = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Arrival orbit"):
            data_arrival = dataset.list_datasets[i].copy()
    x_dep_0 = data_departure[axis_0,:]
    y_dep_1 = data_departure[axis_1,:]
    x_arr_0 = data_arrival[axis_0,:]
    y_arr_1 = data_arrival[axis_1,:]
    
    # Denormalise
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        x_dep_0 *= LU
        y_dep_1 *= LU
        x_arr_0 *= LU
        y_arr_1 *= LU
    
    # Plot orbits
    ax.plot(x_dep_0, y_dep_1,
            label='Departure orbit',
            color=list_colors[0], linestyle=list_linestyles[0])
    ax.plot(x_arr_0, y_arr_1, 
            label='Arrival orbit',
            color=list_colors[1], linestyle=list_linestyles[1])

"""
    Plots a 2D plot from a given dataset.

"""
def plot_2d(dataset):
    
    # Settings
    
    dpi = 200
    
    # Axes
    axis_0 = 0
    axis_1 = 1
    
    # Thrust
    thrust_color = "red"
    
    # System points
    if dataset.dynamical_system.startswith("CR3BP"):
        list_colors_system_points = ["black", "black", "black", "black"]
        list_markers_system_points = ["o", "s", "o", "o"]
        list_plots_system_points = [False, True, True, True]
        
    elif dataset.dynamical_system.startswith("TBP"):
        list_colors_system_points = ["black"]
        list_markers_system_points = ["o"]
        list_plots_system_points = [True]
    
    # Departure arrival
    list_colors_departure_arrival = ["black", "black"]
    list_markers_departure_arrival = ["^", "v"]
    
    # Reference orbits
    list_colors_reference = ["grey", "grey"]
    list_linestyles_references = ["dashed", "dotted"]
    
    # Trajectory
    color_trajectory = "black"
    
    # Normalisation
    denormalise = False

    # Interpolation
    interpolation = True
    interpolation_rate = 15

    # Legend
    show_legend = False
    legend_loc = "lower right"
    
    # Output
    show_grid = True
    save_figure = True
    saving_format = "pdf"
    show_plot = False

    # Get data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i].copy()
            list_names_state = dataset.list_dataset_names[i].copy()
    coord_0 = data_state[axis_0,:]
    coord_1 = data_state[axis_1,:]
    
    # Normalisation
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        coord_0 *= LU
        coord_1 *= LU
        
        # Labels
        list_names_state[axis_0 + 1] = list_names_state[axis_0 + 1].replace(
            "LU", "km")
        list_names_state[axis_1 + 1] = list_names_state[axis_1 + 1].replace(
            "LU", "km")

    if interpolation:
        t_old = np.linspace(0, 1, len(coord_0))  
        t_new = np.linspace(0, 1, interpolation_rate*len(coord_0))  
        coord_0 = interpolate.interp1d(t_old, coord_0, kind='cubic')(t_new)
        coord_1 = interpolate.interp1d(t_old, coord_1, kind='cubic')(t_new)

    # Create plot
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot()
    ax.set_aspect("equal")

    # Set labels
    ax.set_xlabel(list_names_state[axis_0 + 1])
    ax.set_ylabel(list_names_state[axis_1 + 1])
    
    # Plot reference orbits
    plot_reference_orbits(dataset, axis_0, axis_1, ax,
                          list_colors_reference,
                          list_linestyles_references,
                          denormalise)
    
    # Plot departure and arrival points
    plot_departure_arrival(dataset, axis_0, axis_1, ax,
                           list_colors_departure_arrival,
                           list_markers_departure_arrival,
                           denormalise)
    
    # Plot system points
    plot_system_points(dataset, axis_0, axis_1, ax,
                       list_colors_system_points,
                       list_markers_system_points,
                       list_plots_system_points,
                       denormalise)
    
    # Plot Thrust 
    plot_thrust_vector(dataset, axis_0, axis_1,
                       ax, thrust_color,
                       denormalise)
    
    # Plot trajectory
    ax.plot(coord_0, coord_1,
            color=color_trajectory,
            label='Trajectory')
    fig.tight_layout(pad=0.2)
    
    if show_grid:
        plt.grid()
    
    if show_legend: 
        plt.legend(loc=legend_loc)
    
    if save_figure:
        # Get adress
        file_name = dataset.file_name.replace(
            "datasets", "plots")
                
        # Add signature
        signature = ("_2d_" + str(axis_0) + "_" + str(axis_1)
            + "." + saving_format)
        file_name = file_name.replace(
            ".dat", signature)
        
        # Save
        plt.savefig(file_name, bbox_inches='tight')    
       
    if show_plot:   
        plt.show()
    else:
        plt.close(fig)
    
