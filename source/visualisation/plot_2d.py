from classes import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

def plot_system_points(dataset, axis_0, axis_1, ax):
    # TO DO gerer couleurs, taille, etc
    if dataset.dynamical_system.startswith("TBP"):
        ax.scatter(0, 0, label="Central body")
    elif dataset.dynamical_system.startswith("CR3BP"):
        coord_P1 = [-dataset.spacecraft_parameters.constants.MU, 0, 0]
        coord_P2 = [1-dataset.spacecraft_parameters.constants.MU, 0, 0]
        #ax.scatter(coord_P1[axis_0], coord_P1[axis_1], label="$P_1$")
        ax.scatter(coord_P2[axis_0], coord_P2[axis_1], label="$P_2$")
        
    # Departure and Arrival
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i]
    ax.scatter(data_state[axis_0, 0], data_state[axis_1, 0],
               label="Departure")
    ax.scatter(data_state[axis_0, -1], data_state[axis_1, -1],
               label="Arrival")
    
def plot_thrust_vector(dataset, axis_0, axis_1, ax,
                       thrust_scale, thrust_color):
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i]
        elif (dataset.list_dataset_names[i][0] == "Control"):
            data_control = dataset.list_datasets[i]
    coord_0 = data_state[axis_0,:]
    coord_1 = data_state[axis_1,:]
    ucoord_0 = data_control[axis_0,:]
    ucoord_1 = data_control[axis_1,:]
    
    # PLot arrows
    ax.quiver(coord_0[:-1], coord_1[:-1],
              thrust_scale*ucoord_0, thrust_scale*ucoord_1,
              color=thrust_color, label='Thrust')
    
def plot_reference_orbits(dataset, axis_0, axis_1, ax, list_color):
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Departure orbit"):
            data_departure = dataset.list_datasets[i]
        elif (dataset.list_dataset_names[i][0] == "Arrival orbit"):
            data_arrival = dataset.list_datasets[i]
    x_dep_0 = data_departure[axis_0,:]
    y_dep_1 = data_departure[axis_1,:]
    x_arr_0 = data_arrival[axis_0,:]
    y_arr_1 = data_arrival[axis_1,:]
    
    # PLot arrows
    departure_color = list_color[0]
    ax.plot(x_dep_0, y_dep_1,
              color=departure_color, label='Departure orbit')
    arrival_color = list_color[1]
    ax.plot(x_arr_0, y_arr_1, 
            color=arrival_color, label='Arrival orbit')
       

def plot_2d(dataset, axis_0, axis_1):

    # Get data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i]
            list_names_state = dataset.list_dataset_names[i]
    coord_0 = data_state[axis_0,:]
    coord_1 = data_state[axis_1,:]

    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot()
    plt.grid()

    
    # Set labels
    ax.set_xlabel(list_names_state[axis_0 + 1])
    ax.set_ylabel(list_names_state[axis_1 + 1])
    
    # Plot Thrust 
    thrust_scale = 0.3 # TO DO move
    thrust_color = "red" # TO DO move
    plot_thrust_vector(dataset, axis_0, axis_1,
                       ax, thrust_scale, thrust_color)
    
    # Plot reference orbits
    list_colors_reference = ["blue", "orange"]
    plot_reference_orbits(dataset, axis_0, axis_1,
                          ax, list_colors_reference)
    
    # Plot system points
    plot_system_points(dataset, axis_0, axis_1, ax)
    
    # Plot trajectory
    ax.plot(coord_0, coord_1, label='Trajectory')
    
    ax.set_aspect("equal")
    
    plt.legend()
    
    # Show the plot
    plt.show()
    