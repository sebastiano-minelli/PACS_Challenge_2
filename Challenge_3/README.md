### Third challenge  
This project implements the code for the third challenge of PACS course, there are four subfolders:  
- **doc**: it contains just the challenge pdf
- **lib**: it contains the header files
- **src**: it contains *main.cpp* and the *Makefile* to run the program plus the data.txt files
- **files**: it contains the *.vtk* files to be read by paraview  

**How to run**  
The program relies on some PACS course headers. If you have set the environmental variable *PACS_ROOT* to the directory *pacs-example* everything should work fine.  
If something doesn't work then set the variable *PACS_EXAMPLE_PATH* inside *Makefile* to point to your *pacs-example* directory  

This is of course a parallel code, to run it you need to set in the command line the *OpenMP* number of threads and the 
*MPI* processes: *OMP_NUM_THREADS=x mpiexec -np xx ./main*

**Code features**  
- Possibility to run the code in parallel
- Possibility to change BCs and sources by changing the data files
- Possibility to use non-homogeneous BCs
- The output is saved in a *vtk* format and can be read by *Paraview*
- Some tests are performed by running main
  
  
**Files description**
- *Utilities.hpp*: it contains just the dimensions
- *writeVTK.hpp*: it was the header saw during the exercise session changed to use the same containers I used
- *Parameters.hpp*, *ParameterHandler*: they use *muParserXInterface.hpp* to store data, they can be printed out using the suitable function in *ParameterHandler.hpp*
- *InitializeProble.hpp*: it implements the BCs
- *LocalSolver.hpp*: it implements a single iteration of Jacobi
- *JacobiSolver.hpp*: it implements the parallel version of Jacobi
- *main.cpp*: it implements some tests changing the dimension of the matrix, it prints the L2 norms and time

**Useful notes**:  
To perform parallel and sequential scenarios change the parameters on the shell, unfortunately I didn't succed in implementing a script.
The timings seem quite large, probably I'm performing the test wrong
I put two *vtk* files of two pretested problems to ease the visualizations of some features