"""
	main.py

	Purpose: Script to run for visualisation.

	@author Thomas Caleb

	@version 1.0 16/01/2024
    
"""

import os, sys

from classes import Constants, SpacecraftParameters, Dataset
from plot_3d import plot_3d
from plot_2d import plot_2d
from plot_thrust_profile import plot_thrust_profile
from plot_keplerians import plot_keplerian, plot_keplerians
from plot_double_integrator import plot_double_integrator



"""
    Returns a dataset from a filename.

"""
def get_dataset(file_name):
    
    dataset = Dataset()
    dataset.read(file_name)
    return dataset


"""
    Function that runs first.

"""
if __name__ == "__main__":
    # Change directory to DADDy (for windows)
    #os.chdir("../../")
    
    # Get arguments
    list_arguments = sys.argv[1:]
    test_case = int(list_arguments[0])
    T2m_ratio = list_arguments[1]
    ToF = int(list_arguments[2])
    DDP_type = int(list_arguments[3])
    
    extra = ""
    
    file_name = ""
    if test_case == 0:
        file_name = './data/datasets/double_integrator'
    elif test_case == 1:
        file_name = './data/datasets/tbp_SUN_lt_earth_to_mars'
    elif test_case == 2:
        file_name = './data/datasets/tbp_EARTH_lt_leo_to_geo'
    elif test_case == 3:
        file_name = './data/datasets/tbp_EARTH_lt_leo_to_leo'
    elif test_case == 4:
        file_name = './data/datasets/tbp_EARTH_lt_meo_to_meo'
    elif test_case == 5:
        file_name = './data/datasets/cr3bp_EARTH_MOON_lt_haloL2_to_haloL1'
    elif test_case == 6:
        file_name = './data/datasets/cr3bp_EARTH_MOON_lt_nrho_to_dro'
    elif test_case == 7:
        file_name = './data/datasets/cr3bp_EARTH_MOON_lt_lyapunovL1_to_lyapunovL2'
    elif test_case == 8:
        file_name = './data/datasets/cr3bp_EARTH_MOON_lt_dro_to_dro' 
    elif test_case == 9:
        file_name = ''  
    
    # Build name
    if extra != "":
        extra = "_" + extra

    file_name = file_name + "_" + T2m_ratio + "_" + str(ToF) + "_" + str(DDP_type) + extra + ".dat"

    # 
    list_2d = ["tbp", "lyapunov", "dro", "leo", "meo", "halo"]
    list_3d = ["nrho", "halo", "leo", "meo"]
    
    # Load dataset
    dataset = get_dataset(file_name)
    
    # Plots
    if ("double_integrator" in file_name):
        plot_double_integrator(dataset)
    else:   

        plot_thrust_profile(dataset)
        
        for i in list_2d:
            if i in file_name:
                plot_2d(dataset)
                break

        for i in list_3d:
            if i in file_name:
                plot_3d(dataset)
                break

        if ("tbp" in file_name):
            #plot_keplerian(dataset)
            plot_keplerians(dataset)
