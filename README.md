PROGRAMMING ASSIGNMENT PA5
Simulated Annealing
Mark Sears and Navya Sreeram
12/04/2018

==================================================
PROGRAM EXECUTION INSTRUCTIONS

The code assumes that the example design data can be found in a directory "IBM_Benchmarks" (this is set in source code by string "file_path").

Compiling the code can be done with the included makefile. 

To run the code, type in the terminal: "./anneal_std_cell.out param1 param2" 

There are two optional parameters:

param1 is the ibm benchmark to to be used. If no benchmark is specified, the default is "01". The code will look for the design in "./IBM_Benchmarks/ibm##"

param2 is the effort to be used. Effort defines how hard the algorithm will work to find a good solution. Higher efforts lead to lower final wirelength, but longer execution times. If no effort is defined, the default is 6.

==================================================
EXAMPLES

In terminal, navigate to a working directory which includes the source code "anneal_std_cell.cpp", the makefile, and a directory "IBM_Benchmarks".

Type "make" to compile.

Type "./anneal_std_cell.out 01" to run the annealing on benchmark 01.
Type "./anneal_std_cell.out 03 4" to run the annealing on benchmark 03 with effort 4.
Type "./anneal_std_cell.out 12 9" to run the annealing on benchmark 12 with effort 9.

==================================================
DESCRIPTION OF ALGORITHM

The Simulated Annealing method of placement is implemented. Classes for cells, nets and pads are used to organize the design information. The net file is read to store cells, nets and pads information to memory. 
As they are processed, different data structures are created to better organize the information for quick reference and manipulation.

	1) All cells, nets and pads are stored in a vector. This is the primary method of accessing all cells in an arbitrary order.	
	2) A list of cell pointers and pad pointers are used for easy access and computation of wirelengths. This allows the algorithm to compute the wirelengths in minimum amount of time.

The general program flow can be seen in the main function:

	1) The net file is read.
	2) Initial placement of cells and pads are done such that the aspect ratio is approximately 1. The total wirelength is calculated.
	3) A "dry run" of 10,000 cell swaps is performed to obtain the the change in wirelengths and to calculate the standard deviation.
	4) The initial temperature for the algorithm is set proportional to this standard deviation.
	5) The algorithm is executed till the temperature drops to a very low value.

The initial temperature is chosen sufficiently high enough such that a bad move that is 3 times the standard deviation is accepted.

Each step consists of swapping the position of two cells and calculating the change in the wirelength due to the swap. If the new wirelength is smaller, the change is automatically accepted. If the new wirelength is larger, then the probability of accepting is based on the change in wirelength with respect to the temperature.

The algorithm continues performing until the temperature is very low.

==================================================
RESULTS

The netlists range in size from 12,000 cells up to 210,000 cells. More cells means more annealing is required to get similar results. 

Since our implementation scaled the Swaps Per Temperature and initial temperature with the size of the netlist, all benchmarks yielded similar percent change in wirelength. Using an effort of 6, the algorithm achieves a percent change of 82-85%.

