import os
from classes import Constants, SpacecraftParameters, Dataset
from IO import get_dataset
from plot_3d import plot_3d
from plot_2d import plot_2d
from plot_thrust_profile import plot_thrust_profile



if __name__ == "__main__":
    # Change directory to DADDy
    os.chdir("../../")
    
    file_name = './data/datasets/cr3bp_EARTH_MOON_lt_haloL2_to_haloL1.dat'
    #file_name = './data/datasets/tbp_SUN_lt_earth_to_mars.dat'

    dataset = get_dataset(file_name)
    
    # Plot
    plot_3d(dataset)
    plot_2d(dataset)
    plot_thrust_profile(dataset)
