# DADDy
A constrained Differential Algebra-based Differential Dynamic programming solver.

## Requirements
To use the DADDy code, you need a C++11 compiler.

The code was mainly tested on Unix distributions but could be adapted with very low effort.

The only dependency on this code is with the DACE, which can be found at: https://github.com/dacelib/dace.
At installation on the DACE do not forget to use the AlgebraicMatrix<T> class in the CMakeLists.txt.
Moreover, we recommand commenting the 189th line of the file dace/interfaces/cxx/DACEException.cpp too avoid useless and frequent warnings.

## Setting up DADDy
We recommend using CMake to use DADDy and build the test cases.
To use DADDy library, simply clone this repository:
```
git clone https://github.com/ThomasClb/DADDy.git
```
Then create a build directory and run cmake, and then make to compile the examples:
```
cd DADDy
mkdir build
cd build
cmake ../script/test_cases/
make
```
You might need to locate the DACE for cmake if it was not installed using:
```
sudo make install
```

## Running the test case
Nine test cases are given.

In the DADDy folder run:
```
./build/test_case <test_case ID> <parameters> 
```
Leaving the parameters field empty will return an error message that will inform you on what are the available options.

These pararmeters consists in 9 arguments:
	- Test case ID.
	- The adress if the SpacecraftParameter. (They consist in the dynaical system, followed by the thrust to mass ratio of the spacecraft)
	- The DDP method to use. (recommended 2)
	- The number of trajectory segments N?
	- The time of flight in days.
	- Perform fuel optimal optimisation. (1=True, 0=False)
	- Perform solution polishing with Newton method.(1=True, 0=False)
	- Save the data. (1=True, 0=False)
	- Verbosity. (0=Full, 1=Medium, 2=None, 3=Benchmark)

For instance 
	- The double integrator problem from [Lantoine and Russell 2012]:
	```
	./build/test_case 0 ./data/spacecraft_parameters/double_integrator_1e3.dat 2 11 0 0 0 1 0
	```
	- The energy optimal Earth-Mars transfer from [Lantoine and Russell 2012]:
	```
	./build/test_case 1 ./data/spacecraft_parameters/tbp_SUN_5e-4.dat 2 40 348.79 0 1 1 0
	```
	- The fuel optimal Earth-Mars transfer from [Lantoine and Russell 2012]:
	```
	./build/test_case 1 ./data/spacecraft_parameters/tbp_SUN_5e-4.dat 2 40 348.79 1 1 1 0
	```
	- The fuel optimal Halo L2 to Halo L1 transfer from [Aziz et al. 2019]:
	```
	./build/test_case 5 ./data/spacecraft_parameters/cr3bp_EARTH_MOON_5e-4.dat 2 110 20 1 1 1 0
	```
Reading the source code of the test cases is highly recommended to understand how they work and how to use DADDy.

## Visualisation
The visualisation of the results is done using Python 3. Just run:
```
./source/visualization/main.py <test_case ID> <parameters>
```
Where the parameters are:
	- Test case ID.
	- The thrust to mass ratio in N/kg. (string)
	- The time of flight in days. (integer)
	- The DDP method used.
	
For instance 
	- The double integrator problem from [Lantoine and Russell 2012] with DDP method 2:
	```
	./source/visualization/main.py 0 1e3 0 2
	```
	- The fuel optimal Earth-Mars transfer from [Lantoine and Russell 2012] with DDP method 2:
	```
	./source/visualization/main.py 1 5e-4 348 2
	```
	- The fuel optimal Halo L2 to Halo L1 transfer from [Aziz et al. 2019] with DDP method 2:
	```
	./build/test_case 5 5e-4 20 2
	```


