### Second challenge  
This is the folder containing the code for second challenge of PACS course, there are four subfolders:  
- **doc**: it contains just the challenge pdf
- **lib**: it contains the header files
- **src**: it contains *main.cpp* and the *Makefile* to run the program
- **files**: it contains the *.mtx* files to be read

**How to run**  
The program relies on the Chrono.hpp header of the PACS course. If you have set the environmental variable *PACS_ROOT* to the directory *pacs-example* everything should work fine.  
If something doesn't work then set the variable *PACS_EXAMPLE_PATH* inside *Makefile* to point to your *pacs-example* directory  

This is of course a template code, all the variables that can be changed before compiling are declared in the first lines of the *main.cpp* file.  
They are:
- storage ordering method: row-wise or column-wise  
- matrix element type: I set it to double, it should work also with other types (but you should change the read matrices from the files)
- norm types: One, Frobenius or Infinity  

**Code features**  
- Possibility to store matrices row-wise or column-wise (compile time feature)
- Possibility to store matrices in a compress or uncompress format and pass from one state to the other (compressed formats are CSR or CSC depending on the row-wise or column-wise option)
- Matrices can be read from files (by default they are stored in an uncompressed format)
- Matrices can be multiplied by a vector (the vector can be a *std::vector* type or a *Matrix* type)
- There are three possible matrix norms that can be computed: One, Frobenius or Infinity norm
- Methods are programmed to distinguish between the different storing possibilities and implement suitable algorithms based on that

**Useful notes**:  
I wrote some tests in the *main.cpp*.
- Multiplication performance test: it tests if matrix-vector multiplication is faster in the compressed or uncompressed format (one can change the row-wise column-wise storage modifying the variable inside *main.cpp*)
- Multiplication correctness test: it just evaluates matrix-vector multiplication to test if it is correct (the vector in this case is a *Matrix* type)
- Norm test: it just shows all the possible norms for a matrix, actually it computes all the norms both in the compress and uncompress state, the result should be equal in both cases



