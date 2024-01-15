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

def plot_system_points(dataset, ax,
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
            
        if list_plots[0]:
            ax.scatter(coord_center[0], coord_center[1], coord_center[2],
                       label="Central body",
                       color=list_colors[0], marker=list_markers[0])

    elif dataset.dynamical_system.startswith("CR3BP"):
        # Primary coordinates
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
            ax.scatter(coord_P1[0], coord_P1[1], coord_P1[2],
                       label="$P_1$",
                       color=list_colors[0], marker=list_markers[0])
            
        if list_plots[1]:
            ax.scatter(coord_P2[0], coord_P2[1], coord_P2[2],
                       label="$P_2$",
                       color=list_colors[1], marker=list_markers[1])
            
            if list_plots[2]:
                ax.scatter(coord_L1[0], coord_L1[1], coord_L1[2],
                           label="$L_1$",
                           color=list_colors[2], marker=list_markers[2])
                
            if list_plots[3]:
                ax.scatter(coord_L2[0], coord_L2[1], coord_L2[2],
                           label="$L_2$",
                           color=list_colors[3], marker=list_markers[3])

def plot_departure_arrival(dataset, ax,
                           list_colors, list_markers,
                           denormalise):
    # TO DO gerer couleurs, taille, etc
    
    # Retrieve data
    nb_datasets = len(dataset.list_dataset_names)
    for i in range(nb_datasets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i].copy()
            
    # Denormalise
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        data_state *= LU
            
    # Plot departure and arrival
    ax.scatter(data_state[0, 0], data_state[1, 0], data_state[2, 0],
               label="Departure",
               color=list_colors[0], marker=list_markers[0])
    ax.scatter(data_state[0, -1], data_state[1, -1], data_state[2, -1],
               label="Arrival",
               color=list_colors[1], marker=list_markers[1])
    
def plot_thrust_vector(dataset, ax,
                       thrust_scale, thrust_color,
                       denormalise):
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Control"):
            data_control = dataset.list_datasets[i].copy()
    x = data_state[0,:]
    y = data_state[1,:]
    z = data_state[2,:]
    ux = data_control[0,:]
    uy = data_control[1,:]
    uz = data_control[2,:]
    
    # Denormalise
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        THRUSTU = LU*dataset.spacecraft_parameters.constants.THRUSTU
        x *= LU
        y *= LU
        z *= LU
        ux *= THRUSTU
        uy *= THRUSTU
        uz *= THRUSTU
    
    # Plot arrows
    ax.quiver(x[:-1], y[:-1], z[:-1], 
              thrust_scale*ux, thrust_scale*uy, thrust_scale*uz,
              color=thrust_color, label='Thrust')
    
def plot_reference_orbits(dataset, ax,
                          list_colors, list_linestyles,
                          denormalise):
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Departure orbit"):
            data_departure = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Arrival orbit"):
            data_arrival = dataset.list_datasets[i].copy()
    x_dep = data_departure[0,:]
    y_dep = data_departure[1,:]
    z_dep = data_departure[2,:]
    x_arr = data_arrival[0,:]
    y_arr = data_arrival[1,:]
    z_arr = data_arrival[2,:]
    
    # Denormalise
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        x_dep *= LU
        y_dep *= LU
        z_dep *= LU
        x_arr *= LU
        y_arr *= LU
        z_arr *= LU

    # Plot orbits
    ax.plot(x_dep, y_dep, z_dep, 
            label='Departure orbit',
            color=list_colors[0], linestyle=list_linestyles[0])
    ax.plot(x_arr, y_arr, z_arr, 
            label='Arrival orbit',
            color=list_colors[1], linestyle=list_linestyles[1])
    
def plot_3d(dataset):
    
    # Settings
    
    # View 
    azim = 78
    elev = 15
   
    # Thrust
    thrust_scale = 0.2
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
    show_plot = True

    # Get data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i].copy()
            list_names_state = dataset.list_dataset_names[i].copy()
    x = data_state[0,:]
    y = data_state[1,:]
    z = data_state[2,:]
    
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        x *= LU
        y *= LU
        z *= LU
        
        # Labels
        list_names_state[1] = list_names_state[1].replace("LU", "km")
        list_names_state[2] = list_names_state[2].replace("LU", "km")
        list_names_state[3] = list_names_state[3].replace("LU", "km")

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
    plot_thrust_vector(dataset, ax,
                       thrust_scale, thrust_color,
                       denormalise)
    
    # Plot reference orbits
    plot_reference_orbits(dataset, ax,
                          list_colors_reference,
                          list_linestyles_references,
                          denormalise)
    
    # Plot departure and arrival points
    plot_departure_arrival(dataset, ax,
                           list_colors_departure_arrival,
                           list_markers_departure_arrival,
                           denormalise)
    
    # Plot system points
    plot_system_points(dataset, ax,
                       list_colors_system_points,
                       list_markers_system_points,
                       list_plots_system_points,
                       denormalise)
    
    # Plot trajectory
    ax.plot(x, y, z, color=color_trajectory, label='Trajectory')
    
    # View
    ax.view_init(azim=azim, elev=elev)
    
    if show_grid:
        plt.grid()
    
    if show_legend: 
        plt.legend(loc=legend_loc)
    
    if save_figure:
        # Get adress
        file_name = dataset.file_name.replace(
            "datasets", "plots")
                
        # Add signature
        signature = ("_3d" + "." + saving_format)
        file_name = file_name.replace(
            ".dat", signature)
        
        # Save
        plt.savefig(file_name)    
       
    if show_plot:   
        plt.show()
    else:
        plt.close(fig)
    