from classes import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

def plot_system_points(dataset, ax,
                       list_colors, list_markers,
                       list_plots):
    # TO DO gerer couleurs, marker, quel marker plot
    # TO DO point de lagrange
    
    # Discriminate dynamical system
    if dataset.dynamical_system.startswith("TBP"):
        ax.scatter(0, 0, 0,
                   label="Central body")
    elif dataset.dynamical_system.startswith("CR3BP"):
        #ax.scatter(-dataset.spacecraft_parameters.constants.MU, 0, 0,
        #           label="$P_1$")
        ax.scatter(1-dataset.spacecraft_parameters.constants.MU, 0, 0,
                   label="$P_2$")

def plot_plot_departure_arrival(dataset, ax,
                                list_colors, list_markers):
    # TO DO gerer couleurs, taille, etc
    
    # Retrieve data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i]
            
    # Plot departure and arrival
    ax.scatter(data_state[0, 0], data_state[1, 0], data_state[2, 0],
               label="Departure")
    ax.scatter(data_state[0, -1], data_state[1, -1], data_state[2, -1],
               label="Arrival")
    
    
def plot_thrust_vector(dataset, ax,
                       thrust_scale, thrust_color):
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i]
        elif (dataset.list_dataset_names[i][0] == "Control"):
            data_control = dataset.list_datasets[i]
    x = data_state[0,:]
    y = data_state[1,:]
    z = data_state[2,:]
    ux = data_control[0,:]
    uy = data_control[1,:]
    uz = data_control[2,:]
    
    # Plot arrows
    ax.quiver(x[:-1], y[:-1], z[:-1], 
              thrust_scale*ux, thrust_scale*uy, thrust_scale*uz,
              color=thrust_color, label='Thrust')
    
def plot_reference_orbits(dataset, ax,
                          list_colors, list_linestyles):
    # TO DO gerer couleurs, linestyle,faire un template ?

    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Departure orbit"):
            data_departure = dataset.list_datasets[i]
        elif (dataset.list_dataset_names[i][0] == "Arrival orbit"):
            data_arrival = dataset.list_datasets[i]
    x_dep = data_departure[0,:]
    y_dep = data_departure[1,:]
    z_dep = data_departure[2,:]
    x_arr = data_arrival[0,:]
    y_arr = data_arrival[1,:]
    z_arr = data_arrival[2,:]
    departure_color = list_colors[0]
    arrival_color = list_colors[1]
    
    # Plot orbits
    ax.plot(x_dep, y_dep, z_dep, 
              color=departure_color, label='Departure orbit')
    ax.plot(x_arr, y_arr, z_arr, 
            color=arrival_color, label='Arrival orbit')
    

def plot_3d(dataset):

    # Get data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i]
            list_names_state = dataset.list_dataset_names[i]
    x = data_state[0,:]
    y = data_state[1,:]
    z = data_state[2,:]

    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Create cubic bounding box to simulate equal aspect ratio
    max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(x.max()+x.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(y.max()+y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(z.max()+z.min())
    for xb, yb, zb in zip(Xb, Yb, Zb):
       ax.plot([xb], [yb], [zb], 'w')
       
    # White background
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    
    # Set labels
    ax.set_xlabel(list_names_state[1])
    ax.set_ylabel(list_names_state[2])
    ax.set_zlabel(list_names_state[3])
    
    # Plot Thrust 
    thrust_scale = 0.3 # TO DO move
    thrust_color = "red" # TO DO move
    plot_thrust_vector(dataset, ax,
                       thrust_scale, thrust_color)
    
    # Plot reference orbits
    list_colors_reference = ["blue", "orange"]
    list_linestyles_references = ["dotted", "dotted"]
    plot_reference_orbits(dataset, ax,
                          list_colors_reference,
                          list_linestyles_references)
    
    # Plot system points
    list_colors_system_points = ["black"]
    list_markers_system_points = ["o"]
    list_plots_system_points = ["Center"]
    plot_system_points(dataset, ax,
                       list_colors_system_points,
                       list_markers_system_points,
                       list_plots_system_points)
    
    # Plot trajectory
    color_trajectory = "green"
    ax.plot(x, y, z, color=color_trajectory, label='Trajectory')
    
    # Save
    # TO DO
    
    # Show the plot
    plt.show()
    