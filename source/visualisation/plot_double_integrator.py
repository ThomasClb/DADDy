"""
	plot_double_integrator.py

	Purpose: Implements the functions to produce plot for double_integrator.

	@author Thomas Caleb

	@version 1.0 26/04/2024
    
"""

import numpy as np
import matplotlib.pyplot as plt
#plt.switch_backend('agg')
import scipy.interpolate as interpolate

from misc import get_Lagrange_point
from classes import Dataset



"""
    Plots a from a given dataset.

"""
def plot_double_integrator_u(dataset):
    
    # Settings
    dpi = 200
        
    # Departure arrival
    list_colors = ["blue", "green", "purple"]
    list_markers = ["o", "s", "^"]

    # Legend
    show_legend = True
    legend_loc_control = "upper left"
    legend_loc_state = "lower left"
    
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
        if (dataset.list_dataset_names[i][0] == "Control"):
            data_control = dataset.list_datasets[i].copy()
            list_names_control = dataset.list_dataset_names[i].copy()
    coord_0 = data_state[0,:]
    coord_1 = data_state[1,:]
    coord_2 = data_state[2,:]
    control_0 = data_control[0,:]
    control_1 = data_control[1,:]
    control_2 = data_control[2,:]
    N = len(control_1)

    # Change labels
    list_names_state[0 + 1] = list_names_state[0 + 1].replace(
        "[LU]", "")
    list_names_state[1 + 1] = list_names_state[1 + 1].replace(
        "[LU]", "")
    list_names_state[2 + 1] = list_names_state[2 + 1].replace(
        "[LU]", "")
    list_names_control[0 + 1] = list_names_control[0 + 1].replace(
        "[THRUSTU]", "")
    list_names_control[1 + 1] = list_names_control[1 + 1].replace(
        "[THRUSTU]", "")
    list_names_control[2 + 1] = list_names_control[2 + 1].replace(
        "[THRUSTU]", "")


    # Create plot
    fig, (ax1) = plt.subplots(1, 1, dpi=dpi)

    # Set labels
    ax1.set_xlabel("Stage [-]")
    ax1.set_ylabel("Controls [-]")
    
    # Plot control
    ax1.scatter(range(N), control_0,
            color=list_colors[0],
            marker=list_markers[0],
            label=list_names_control[1])
    ax1.scatter(range(N), control_1,
            color=list_colors[1],
            marker=list_markers[1],
            label=list_names_control[2][1:])
    ax1.scatter(range(N), control_2,
            color=list_colors[2],
            marker=list_markers[2],
            label=list_names_control[3][1:])
    fig.tight_layout(pad=0.2)
    
    if show_grid:
        ax1.grid()
    
    if show_legend: 
        ax1.legend(loc=legend_loc_control)
    
    if save_figure:
        # Get adress
        file_name = dataset.file_name.replace(
            "datasets", "plots")
                
        # Add signature
        signature = ("_u." + saving_format)
        file_name = file_name.replace(
            ".dat", signature)
        
        # Save
        plt.savefig(file_name, bbox_inches='tight')    
       
    if show_plot:   
        plt.show()
    else:
        plt.close(fig)
        
def plot_double_integrator_x(dataset):
    
    # Settings
    dpi = 200
        
    # Departure arrival
    list_colors = ["blue", "green", "purple"]
    list_markers = ["o", "s", "^"]

    # Legend
    show_legend = True
    legend_loc_control = "upper left"
    legend_loc_state = "lower left"
    
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
        if (dataset.list_dataset_names[i][0] == "Control"):
            data_control = dataset.list_datasets[i].copy()
            list_names_control = dataset.list_dataset_names[i].copy()
    coord_0 = data_state[0,:]
    coord_1 = data_state[1,:]
    coord_2 = data_state[2,:]
    control_0 = data_control[0,:]
    control_1 = data_control[1,:]
    control_2 = data_control[2,:]
    N = len(control_1)

    # Change labels
    list_names_state[0 + 1] = list_names_state[0 + 1].replace(
        "[LU]", "")
    list_names_state[1 + 1] = list_names_state[1 + 1].replace(
        "[LU]", "")
    list_names_state[2 + 1] = list_names_state[2 + 1].replace(
        "[LU]", "")
    list_names_control[0 + 1] = list_names_control[0 + 1].replace(
        "[THRUSTU]", "")
    list_names_control[1 + 1] = list_names_control[1 + 1].replace(
        "[THRUSTU]", "")
    list_names_control[2 + 1] = list_names_control[2 + 1].replace(
        "[THRUSTU]", "")


    # Create plot
    fig, ax2 = plt.subplots(1, 1, dpi=dpi)

    # Set labels
    ax2.set_xlabel("Stage [-]")
    ax2.set_ylabel("States [-]")
    
    # Plot control

    # Plot control
    ax2.scatter(range(N+1), coord_0,
            color=list_colors[0],
            marker=list_markers[0],
            label=list_names_state[1])
    ax2.scatter(range(N+1), coord_1,
            color=list_colors[1],
            marker=list_markers[1],
            label=list_names_state[2][1:])
    ax2.scatter(range(N+1), coord_2,
            color=list_colors[2],
            marker=list_markers[2],
            label=list_names_state[3][1:])
    fig.tight_layout(pad=0.2)

    if show_grid:
        ax2.grid()
    
    if show_legend: 
        ax2.legend(loc=legend_loc_state)
    
    if save_figure:
        # Get adress
        file_name = dataset.file_name.replace(
            "datasets", "plots")
                
        # Add signature
        signature = ("_x." + saving_format)
        file_name = file_name.replace(
            ".dat", signature)
        
        # Save
        plt.savefig(file_name, bbox_inches='tight')    
       
    if show_plot:   
        plt.show()
    else:
        plt.close(fig)
    


def plot_double_integrator(dataset):
    
    plot_double_integrator_x(dataset)
    plot_double_integrator_u(dataset)
    
    
    # Settings
    dpi = 200
        
    # Departure arrival
    list_colors = ["blue", "green", "purple"]
    list_markers = ["o", "s", "^"]

    # Legend
    show_legend = True
    legend_loc_control = "upper left"
    legend_loc_state = "lower left"
    
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
        if (dataset.list_dataset_names[i][0] == "Control"):
            data_control = dataset.list_datasets[i].copy()
            list_names_control = dataset.list_dataset_names[i].copy()
    coord_0 = data_state[0,:]
    coord_1 = data_state[1,:]
    coord_2 = data_state[2,:]
    control_0 = data_control[0,:]
    control_1 = data_control[1,:]
    control_2 = data_control[2,:]
    N = len(control_1)

    # Change labels
    list_names_state[0 + 1] = list_names_state[0 + 1].replace(
        "[LU]", "")
    list_names_state[1 + 1] = list_names_state[1 + 1].replace(
        "[LU]", "")
    list_names_state[2 + 1] = list_names_state[2 + 1].replace(
        "[LU]", "")
    list_names_control[0 + 1] = list_names_control[0 + 1].replace(
        "[THRUSTU]", "")
    list_names_control[1 + 1] = list_names_control[1 + 1].replace(
        "[THRUSTU]", "")
    list_names_control[2 + 1] = list_names_control[2 + 1].replace(
        "[THRUSTU]", "")


    # Create plot
    fig, (ax1, ax2) = plt.subplots(1, 2, dpi=dpi)

    # Set labels
    ax1.set_xlabel("Stage [-]")
    ax1.set_ylabel("Controls [-]")
    ax2.set_xlabel("Stage [-]")
    ax2.set_ylabel("States [-]")
    
    # Plot control
    ax1.scatter(range(N), control_0,
            color=list_colors[0],
            marker=list_markers[0],
            label=list_names_control[1])
    ax1.scatter(range(N), control_1,
            color=list_colors[1],
            marker=list_markers[1],
            label=list_names_control[2][1:])
    ax1.scatter(range(N), control_2,
            color=list_colors[2],
            marker=list_markers[2],
            label=list_names_control[3][1:])

    # Plot control
    ax2.scatter(range(N+1), coord_0,
            color=list_colors[0],
            marker=list_markers[0],
            label=list_names_state[1])
    ax2.scatter(range(N+1), coord_1,
            color=list_colors[1],
            marker=list_markers[1],
            label=list_names_state[2][1:])
    ax2.scatter(range(N+1), coord_2,
            color=list_colors[2],
            marker=list_markers[2],
            label=list_names_state[3][1:])
    fig.tight_layout(pad=0.2)
    
    if show_grid:
        ax1.grid()
        ax2.grid()
    
    if show_legend: 
        ax1.legend(loc=legend_loc_control)
        ax2.legend(loc=legend_loc_state)
    
    if save_figure:
        # Get adress
        file_name = dataset.file_name.replace(
            "datasets", "plots")
                
        # Add signature
        signature = ("." + saving_format)
        file_name = file_name.replace(
            ".dat", signature)
        
        # Save
        plt.savefig(file_name, bbox_inches='tight')    
       
    if show_plot:   
        plt.show()
    else:
        plt.close(fig)
    
