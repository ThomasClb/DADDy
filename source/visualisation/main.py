from classes import Constants, SpacecraftParameters, Dataset
from IO import get_dataset
from plot_3d import plot_3d
from plot_2d import plot_2d
from plot_thrust_profile import plot_thrust_profile



if __name__ == "__main__":
    file_name = '../../data/datasets/cr3bp_halo_to_halo.dat'
    dataset = get_dataset(file_name)
    
    # Plot
    #plot_3d(dataset)
    plot_2d(dataset, 1, 2)
    
    #plot_thrust_profile(dataset)
