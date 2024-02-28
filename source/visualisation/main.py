"""
	main.py

	Purpose: Script to run for visualisation.

	@author Thomas Caleb

	@version 1.0 16/01/2024
    
"""

import os

from classes import Constants, SpacecraftParameters, Dataset
from plot_3d import plot_3d
from plot_2d import plot_2d
from plot_thrust_profile import plot_thrust_profile
from plot_keplerians import plot_keplerian, plot_keplerians



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
    # Change directory to DADDy
    os.chdir("../../")
    
    
    DDP_type = "_2"
    T2m_ratio = "_5e-4"
        
    #file_name = './data/datasets/tbp_SUN_lt_earth_to_mars'
    #file_name = './data/datasets/tbp_EARTH_lt_leo_to_geo'
    #file_name = './data/datasets/tbp_EARTH_lt_leo_to_leo'
    file_name = './data/datasets/tbp_EARTH_lt_meo_to_meo'
    #file_name = './data/datasets/cr3bp_EARTH_MOON_lt_haloL2_to_haloL1'
    #file_name = './data/datasets/cr3bp_EARTH_MOON_lt_nrho_to_dro'
    #file_name = './data/datasets/cr3bp_EARTH_MOON_lt_lyapunovL1_to_lyapunovL2'
    #file_name = './data/datasets/cr3bp_EARTH_MOON_lt_dro_to_dro'
    
    
    
    file_name = file_name + DDP_type + T2m_ratio + ".dat"

    
    # Load dataset
    dataset = get_dataset(file_name)
    
    # Plots
    plot_3d(dataset)
    plot_2d(dataset)
    plot_thrust_profile(dataset)
    #plot_keplerian(dataset)
    plot_keplerians(dataset)
