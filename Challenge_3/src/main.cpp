#include "ParameterHandler.hpp"
#include "writeVTK.hpp"
#include <iostream>
#include <vector>
#include "mpi_utils.hpp"
#include "SafeMPI.hpp"
#include "JacobiSolver.hpp"

int main(int argc, char **argv)
{
    param::ParameterHandler params("data.txt");
    params.show_data();

    const int N = params.coefficients.n;

    const double h = 1.0/ N;

    std::vector<double> exacSol(N * N, 0.0);
    JacobiSolver jac_solver(exacSol, params, argc, argv);
    MPI_Init(&argc, &argv);
    auto [sol, norm, n_it, rank] = jac_solver.solve();
    
    if(rank == 0)
    {
        std::cout << "Local norm: " << norm << std::endl;
        std::cout << "Number of iterations: " << n_it << std::endl;

        generateVTKFile("../files/output.vtk", sol, N, N, h, h);
    }
    MPI_Finalize();
    


    return 0;
}