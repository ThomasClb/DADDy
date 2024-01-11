from classes import Dataset
import numpy as np
import matplotlib.pyplot as plt

def plot_stairs_profile(dt, u, ax, label, color):
    nb_point_stair = 100
    nb_point = len(dt)
    t = []
    u_stairs = []
    current_time = 0
    for i in range(nb_point):
       
        dt_i = dt[i]
        u_i = u[i]
        for j in range(nb_point_stair):
            t.append(current_time)
            u_stairs.append(u_i)
            current_time = current_time + dt_i*j/(nb_point_stair - 1)
            
    ax.plot(t, u_stairs, label=label, color=color)
    
def plot_thrust_profile(dataset):
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i]
        elif (dataset.list_dataset_names[i][0] == "Control"):
            data_control = dataset.list_datasets[i]
    ux = data_control[0,:]
    uy = data_control[1,:]
    uz = data_control[2,:]
    dt = data_state[7,:-1]
    u = np.sqrt(ux*ux + uy*uy + uz*uz)
    max_trust = dataset.spacecraft_parameters.thrust

    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot()
    plt.grid()

    # Set labels
    ax.set_xlabel("Time [TU]")
    ax.set_ylabel("Thrust norm [THRUSTU]")
    
    # Plot Thrust 
    color = "blue" # TO DO move
    label = "Thrust"
    plot_stairs_profile(dt, u, ax, label, color)
    
    plot_components = False
    if plot_components:
        color = "green" # TO DO move
        label = "ux"
        plot_stairs_profile(dt, ux, ax, label, color)
        
        color = "cyan" # TO DO move
        label = "uy"
        plot_stairs_profile(dt, uy, ax, label, color)
        
        color = "grey" # TO DO move
        label = "uz"
        plot_stairs_profile(dt, uz, ax, label, color)
    
        
    # Plot Thrust 
    color = "red" # TO DO move
    label = "Max thrust"
    plot_stairs_profile(dt, u*0 + max_trust, ax, label, color)
    
    # Show legend
    plt.legend()
    
    # Show the plot
    plt.show()
    