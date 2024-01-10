from classes import Constants, SpacecraftParameters, Dataset

def get_dataset(file_name):
    
    # Read constants and spacecraft parameters
    dataset = Dataset()
    dataset.read(file_name)
    
    
    
    # Check     
    print("Attribute 1:", dataset.dynamical_system)
    print("Attribute 2:", dataset.spacecraft_parameters.constants.MU)
    print("Attribute 3:", dataset.spacecraft_parameters.constants.WU)
    print("Attribute 4:", dataset.spacecraft_parameters.constants.MASSU)
    print("Attribute 5:", dataset.spacecraft_parameters.constants.TU)
    print("Attribute 6:", dataset.spacecraft_parameters.constants.VU)
    print("Attribute 7:", dataset.spacecraft_parameters.constants.THRUSTU)
    print("Attribute 8:", dataset.spacecraft_parameters.initial_mass)
    print("Attribute 9:", dataset.spacecraft_parameters.dry_mass)
    print("Attribute 10:", dataset.spacecraft_parameters.thrust)
    print("Attribute 11:", dataset.spacecraft_parameters.Isp)
    print(dataset.list_dataset_names)
    print(dataset.list_datasets)
    
    return dataset

if __name__ == "__main__":
    file_name = '../../data/datasets/test.dat'
    get_dataset(file_name)
