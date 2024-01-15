import numpy as np

class Constants:
        
    def __init__(self, mu=1, lu=1, wu=1, massu=1000):
        self.MU = float(mu)
        self.LU = float(lu)
        self.WU = float(wu)
        self.MASSU = float(massu)
        self.TU = 1/self.WU
        self.VU = self.LU*self.WU
        self.THRUSTU = 1000 * self.VU * self.MASSU * self.WU;
        
    def read(self , file):
        file.readline() # skip first
        self.MU = float(file.readline().strip())
        self.LU = float(file.readline().strip())
        self.WU = float(file.readline().strip())
        self.MASSU = float(file.readline().strip())
        self.TU = 1/self.WU
        self.VU = self.LU*self.WU
        self.THRUSTU = 1000 * self.VU * self.MASSU * self.WU;

class SpacecraftParameters:
    
    def __init__(self, constants=Constants(),
                 initial_mass=1, dry_mass=0.5, thrust=1, Isp=1):
        self.constants = constants
        self.initial_mass = float(initial_mass)
        self.dry_mass = float(dry_mass)
        self.thrust = float(thrust)
        self.Isp = float(Isp)

    def read(self , file):
        self.constants.read(file)
        file.readline() # skip first
        self.initial_mass = float(file.readline().strip())
        self.dry_mass = float(file.readline().strip())
        self.thrust = float(file.readline().strip())
        self.Isp = float(file.readline().strip())
        
class Dataset:
    
    def __init__(self,
                 file_name="",
                 dynamical_system="",
                 spacecraft_parameters=SpacecraftParameters(),
                 list_dataset_names=[],
                 list_datasets=[]):
        self.file_name = file_name
        self.dynamical_system = dynamical_system
        self.spacecraft_parameters = spacecraft_parameters
        self.list_dataset_names = list_dataset_names
        self.list_datasets = list_datasets

    def read(self , file_name):
        with open(file_name, 'r') as file:
        
            # Read header
            self.file_name = file.readline().strip()
            self.dynamical_system = file.readline().strip()
            self.spacecraft_parameters.read(file)
            
            # Get number of datasets
            nb_datasets = int(file.readline().strip())
            
            list_dataset_names = []
            list_datasets = []
            for i in range(nb_datasets):
                # Init dataset names
                list_dataset_names_i = []
                list_dataset_names_i.append(file.readline().strip())
                
                # Get array size, rows and colomn are inverted
                size_array = file.readline().strip().split(",")
                nb_rows = int(size_array[1]) 
                nb_cols = int(size_array[0])
                dataset = np.zeros((nb_rows, nb_cols))
                
                # Get names
                list_dataset_names_i = list_dataset_names_i + file.readline().strip().split(",")
                list_dataset_names_i.pop(nb_rows + 1)
                list_dataset_names.append(list_dataset_names_i)
                
                # Get data
                for k in range(nb_cols):
                    # Get data
                    list_data = file.readline().strip().split(",")
                    for l in range(nb_rows):
                        dataset[l, k] = float(list_data[l])
                list_datasets.append(dataset)
            self.list_dataset_names = list_dataset_names;
            self.list_datasets = list_datasets;
                        