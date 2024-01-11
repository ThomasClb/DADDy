from classes import Constants, SpacecraftParameters, Dataset

def get_dataset(file_name):
    
    # Read constants and spacecraft parameters
    dataset = Dataset()
    dataset.read(file_name)
    return dataset
