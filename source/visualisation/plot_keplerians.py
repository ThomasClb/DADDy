"""
	plot_keplerians.py

	Purpose: Implements the functions to output Keplerian visualisations

	@author Thomas Caleb

	@version 1.0 25/01/2024
    
"""

import numpy as np
import matplotlib.pyplot as plt
#plt.switch_backend('agg')
from classes import Dataset

"""convert cartesian to orbital parameters"""
def cart_2_kep(x, y, z, vx, vy, vz, mu):
    # Compute state, velocity and angular momentum
    r_vector = np.array([x, y, z])
    v_vector = np.array([vx, vy, vz])
    r = np.linalg.norm(r_vector)
    v = np.linalg.norm(v_vector)
    h_vector = np.cross(r_vector, v_vector)
    h = np.linalg.norm(h_vector)
        
    # Get a, e, inc
    sma = 1/(2/r - v**2/mu)
    ecc_vector = np.cross(v_vector, h_vector)/mu - r_vector/r
    ecc = np.linalg.norm(ecc_vector)
    inc = np.arccos(h_vector[2]/h)
    
    # Get RAAN
    N_vector = np.cross([0, 0, 1], h_vector)
    N = np.linalg.norm(N_vector)
    if N_vector[1] >= 0:
        if N != 0:
            RAAN = np.arccos(N_vector[0]/N)
        else:
            RAAN = np.pi/2
    else:
        if N != 0:
            RAAN = 2*np.pi - np.arccos(N_vector[0]/N)
        else:
            RAAN = 1.5*np.pi
        
    # Get nu
    if ecc != 0:
        
        # Get omega
        alpha = np.dot(ecc_vector, N_vector)
        beta = ecc*N
        if ecc_vector[2] >= 0:
            if beta != 0:
                omega = np.arccos(max(-1, min(alpha/beta, 1)))
            else:
                omega = np.pi/2
        else:
            if beta != 0:
                omega = 2*np.pi - np.arccos(max(-1, min(alpha/beta, 1)))
            else:
                omega = 1.5*np.pi

        # Get nu
        alpha = np.dot(ecc_vector, r_vector)
        beta = ecc*r
        if np.dot(r_vector, v_vector) >= 0:
            if beta != 0:
                nu = np.arccos(max(-1, min(alpha/beta, 1)))
            else:
                nu = np.pi/2
        else:
            if beta != 0:
                nu = 2*np.pi - np.arccos(max(-1, min(alpha/beta, 1)))
            else:
                nu = 1.5*np.pi

    else:
        omega = 0
        if y >= 0:
            nu = np.arccos(x/r)
        else:
            nu = 2*np.pi - np.arccos(x/r)
            
    return sma, ecc, inc, RAAN, omega, nu


"""
    Plots an orbital element for a given transfer dataset.

"""
def plot_keplerian(dataset):
    # Settings
    
    dpi = 200
    
    # Keplerian selection
    coordinate = 2
    color = "black"
    
    # Normalisation
    denormalise = True
    except_sma = True
    
    # Legend
    show_legend = False
    legend_loc = "lower right"

    # Output
    show_grid = True
    save_figure = True
    saving_format = "pdf"
    show_plot = False
    
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i].copy()
    x = data_state[0,:]
    y = data_state[1,:]
    z = data_state[2,:]
    vx = data_state[3,:]
    vy = data_state[4,:]
    vz = data_state[5,:]
    dt = data_state[7,:]
    mu = dataset.spacecraft_parameters.constants.MU/dataset.spacecraft_parameters.constants.MU
    
    # Convert data to keplerians
    sma = np.zeros(x.shape)
    ecc = np.zeros(x.shape)
    inc = np.zeros(x.shape)
    RAAN = np.zeros(x.shape)
    omega = np.zeros(x.shape)
    nu = np.zeros(x.shape)
    t = np.zeros(x.shape)
    for i in range(len(x)):
        sma[i], ecc[i], inc[i], RAAN[i], omega[i], nu[i] = cart_2_kep(
            x[i], y[i], z[i], vx[i], vy[i], vz[i], mu)
        if i > 0:
            t[i] = t[i - 1] + dt[i]

    # Get label
    x_label = "Time [TU]"
    if coordinate == 0:
        y_label = "Semi-major axis [LU]"
    elif coordinate == 1:
        y_label = "Eccentricity [-]"
    elif coordinate == 2:
        y_label = "Inclination [rad]"
    elif coordinate == 3:
        y_label = "RAAN [rad]"
    elif coordinate == 4:
        y_label = "Argument of periapsis [rad]"
    elif coordinate == 5:
        y_label = "True anomaly [rad]"

    # Normalisation
    rad2deg = 180/np.pi
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        TU = dataset.spacecraft_parameters.constants.TU
        inc *= rad2deg
        RAAN *= rad2deg
        omega *= rad2deg
        nu *= rad2deg
        t *= (TU/86400)
        if not except_sma:
            sma *= LU
        
        # Labels
        x_label = x_label.replace("TU", "days")
        if ("[km]" in y_label and not except_sma):
            y_label = y_label.replace("km", "LU")
        elif ("[rad]" in y_label):
            y_label = y_label.replace("rad", "deg")
    
    # Select plot
    if coordinate == 0:
        y = sma
    elif coordinate == 1:
        y = ecc
    elif coordinate == 2:
        y = inc
    elif coordinate == 3:
        y = RAAN
    elif coordinate == 4:
        y = omega
    elif coordinate == 5:
        y = nu

    # Create plot
    fig, ax = plt.subplots(dpi=dpi)

    # Make plot
    ax.plot(t, y, color=color)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    fig.tight_layout(pad=0.3)

    if show_grid:
        plt.grid()
    
    if show_legend: 
        plt.legend(loc=legend_loc)
    
    if save_figure:
        # Get adress
        file_name = dataset.file_name.replace(
            "datasets", "plots")
                
        # Add signature
        signature = ("_keplerian_" + str(coordinate) + "." + saving_format)
        file_name = file_name.replace(
            ".dat", signature)
        
        # Save
        plt.savefig(file_name, bbox_inches='tight')    
       
    if show_plot:   
        plt.show()
    else:
        plt.close(fig)
        
"""
    Plots a an orbital element for a given transfer dataset.

"""
def plot_keplerians(dataset):
    # Settings
    
    dpi = 200
    
    # Keplerian selection
    list_coordinate = [0, 1, 2]
    color = "black"
    
    # Normalisation
    denormalise = True
    except_sma = True
    
    # Legend
    show_legend = False
    legend_loc = "lower right"

    # Output
    show_grid = True
    save_figure = True
    saving_format = "pdf"
    show_plot = False
    
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "State"):
            data_state = dataset.list_datasets[i].copy()
    x = data_state[0,:]
    y = data_state[1,:]
    z = data_state[2,:]
    vx = data_state[3,:]
    vy = data_state[4,:]
    vz = data_state[5,:]
    dt = data_state[7,:]
    mu = dataset.spacecraft_parameters.constants.MU/dataset.spacecraft_parameters.constants.MU
    
    # Convert data to keplerians
    sma = np.zeros(x.shape)
    ecc = np.zeros(x.shape)
    inc = np.zeros(x.shape)
    RAAN = np.zeros(x.shape)
    omega = np.zeros(x.shape)
    nu = np.zeros(x.shape)
    t = np.zeros(x.shape)
    for i in range(len(x)):
        sma[i], ecc[i], inc[i], RAAN[i], omega[i], nu[i] = cart_2_kep(
            x[i], y[i], z[i], vx[i], vy[i], vz[i], mu)
        if i > 0:
            t[i] = t[i - 1] + dt[i]

    # Get label
    x_label = "Time [TU]"
    list_y_label = []
    if 0 in list_coordinate:
        list_y_label.append("Semi-major axis [LU]")
    if 1 in list_coordinate:
       list_y_label.append("Eccentricity [-]")
    if 2 in list_coordinate:
        list_y_label.append("Inclination [rad]")
    if 3 in list_coordinate:
        list_y_label.append("RAAN [rad]")
    if 4 in list_coordinate:
        list_y_label.append("Argument of periapsis [rad]")
    if 5 in list_coordinate:
        list_y_label.append("True anomaly [rad]")

    # Normalisation
    rad2deg = 180/np.pi
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        TU = dataset.spacecraft_parameters.constants.TU
        inc *= rad2deg
        RAAN *= rad2deg
        omega *= rad2deg
        nu *= rad2deg
        t *= (TU/86400)
        if not except_sma:
            sma *= LU
        
        # Labels
        x_label = x_label.replace("TU", "days")
        for i, y_label in enumerate(list_y_label):
            if ("[km]" in y_label and not except_sma):
                list_y_label[i] = y_label.replace("km", "LU")
            elif ("[rad]" in y_label):
                list_y_label[i] = y_label.replace("rad", "deg")
    
    # Select plot
    list_y = []
    if 0 in list_coordinate:
        list_y.append(sma)
    if 1 in list_coordinate:
        list_y.append(ecc)
    if 2 in list_coordinate:
        list_y.append(inc)
    if 3 in list_coordinate:
        list_y.append(RAAN)
    if 4 in list_coordinate:
        list_y.append(omega)
    if 5 in list_coordinate:
        list_y.append(nu)

    # Create plot TO DO shapes    
    if len(list_coordinate) == 1:
        list_index = [(0)]
        shape = (len(list_coordinate), 1)
        fig, ax = plt.subplots(shape[0], shape[1],dpi=dpi)
        fig.tight_layout(pad=0.3)
        
    elif len(list_coordinate) == 2:
        list_index = [(0), (1)]
        shape = (len(list_coordinate), 1)
        fig, ax = plt.subplots(shape[0], shape[1],dpi=dpi)
        fig.tight_layout(pad=0.3)

    elif len(list_coordinate) == 3:
        list_index = [(0), (1), (2)]
        shape = (len(list_coordinate), 1)
        fig, ax = plt.subplots(shape[0], shape[1],dpi=dpi)
        fig.tight_layout(pad=0.3)
        
    elif len(list_coordinate) == 4:
        list_index = [(0, 0), (1, 0), (0, 1), (1, 1)]
        shape = (2, 2)
        fig, ax = plt.subplots(shape[0], shape[1],dpi=dpi)
        fig.tight_layout(pad=0.3)
        
    elif len(list_coordinate) == 5:
        list_index = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1)]
        shape = (3, 2)
        fig, ax = plt.subplots(shape[0], shape[1],dpi=dpi)
        fig.tight_layout(pad=0.3)
        
    elif len(list_coordinate) == 6:
        list_index = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)]
        shape = (3, 2)
        fig, ax = plt.subplots(shape[0], shape[1],dpi=dpi)
        fig.tight_layout(pad=0.3)

    # Make plot
    for i, y in enumerate(list_y):
        ax_i = ax[list_index[i]]
        ax_i.plot(t, y, color=color)
        ax_i.set_xlabel(x_label)
        ax_i.set_ylabel(list_y_label[i])
    

    if show_grid:
        plt.grid()
    
    if show_legend: 
        plt.legend(loc=legend_loc)
    
    if save_figure:
        # Get adress
        file_name = dataset.file_name.replace(
            "datasets", "plots")
                
        # Add signature
        signature_list = ""
        for coordinate in list_coordinate:
            signature_list = signature_list + "_" + str(coordinate)
        signature = ("_keplerians" + signature_list + "." + saving_format)
        file_name = file_name.replace(
            ".dat", signature)
        
        # Save
        plt.savefig(file_name,bbox_inches='tight')    
       
    if show_plot:   
        plt.show()
    else:
        plt.close(fig)
    