from classes import Dataset
import numpy as np
import matplotlib.pyplot as plt

def get_Lagrange_point(dataset, point):
    q = dataset.spacecraft_parameters.constants.MU
    eps = (q/3.0)**(1/3.0)
    
    if point == 1:  
        x = 1.0-q - eps + eps**2/3 + eps**3/9

    elif point == 2:  
        x = 1.0 - q + eps + eps**2/3 - eps**3/9
    
    return np.array([x, 0.0, 0.0])

def plot_system_points(dataset, axis_0, axis_1, ax,
                       list_colors, list_markers,
                       list_plots,
                       denormalise):    
    # Discriminate dynamical system
    if dataset.dynamical_system.startswith("TBP"):
        # Make points
        coord_center = np.array([0, 0, 0])
        
        # Denormalise
        if denormalise:
            LU = dataset.spacecraft_parameters.constants.LU
            coord_center *= LU
        
        if list_plots[0]:
            ax.scatter(coord_center[axis_0], coord_center[axis_1],
                       label="Central body",
                       color=list_colors[0], marker=list_markers[0])

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
            
        # Plots
        if list_plots[0]:
            ax.scatter(coord_P1[axis_0], coord_P1[axis_1],
                       label="$P_1$",
                       color=list_colors[0], marker=list_markers[0])
        
        if list_plots[1]:
            ax.scatter(coord_P2[axis_0], coord_P2[axis_1],
                       label="$P_2$",
                       color=list_colors[1], marker=list_markers[1])
            
        if list_plots[2]:
            ax.scatter(coord_L1[axis_0], coord_L1[axis_1],
                       label="$L_1$",
                       color=list_colors[2], marker=list_markers[2])
        
        if list_plots[3]:
            ax.scatter(coord_L2[axis_0], coord_L2[axis_1],
                       label="$L_2$",
                       color=list_colors[3], marker=list_markers[3])    
        
        
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
               label="Departure",
               color=list_colors[0], marker=list_markers[0])
    ax.scatter(data_state[axis_0, -1], data_state[axis_1, -1],
               label="Arrival",
               color=list_colors[1], marker=list_markers[1])
    
def plot_thrust_vector(dataset, axis_0, axis_1, ax,
                       thrust_scale, thrust_color,
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
    
    # PLot arrows
    ax.quiver(coord_0[:-1], coord_1[:-1],
              thrust_scale*ucoord_0, thrust_scale*ucoord_1,
              color=thrust_color, label='Thrust')
    
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
       
def plot_2d(dataset):
    
    # Settings
    
    # Axes
    axis_0 = 0
    axis_1 = 1
    
    # Thrust
    thrust_scale = 0.3
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
    
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        coord_0 *= LU
        coord_1 *= LU
        
        # Labels
        list_names_state[axis_0 + 1] = list_names_state[axis_0 + 1].replace(
            "LU", "km")
        list_names_state[axis_1 + 1] = list_names_state[axis_1 + 1].replace(
            "LU", "km")

    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_aspect("equal")

    # Set labels
    ax.set_xlabel(list_names_state[axis_0 + 1])
    ax.set_ylabel(list_names_state[axis_1 + 1])
    
    # Plot Thrust 
    plot_thrust_vector(dataset, axis_0, axis_1,
                       ax, thrust_scale, thrust_color,
                       denormalise)
    
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
    
    # Plot trajectory
    ax.plot(coord_0, coord_1,
            color=color_trajectory,
            label='Trajectory')
    
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
        plt.savefig(file_name)    
       
    if show_plot:   
        plt.show()
    else:
        plt.close(fig)
    