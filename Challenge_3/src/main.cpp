#include "ParameterHandler.hpp"
#include "writeVTK.hpp"
#include <iostream>
#include <vector>
#include "mpi_utils.hpp"
#include "SafeMPI.hpp"
#include "JacobiSolver.hpp"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    param::ParameterHandler params("data.txt");
    if(rank == 0)
        params.show_data();

    const int N = params.coefficients.n;

    const double h = 1.0/ N;

    std::vector<double> exacSol(N * N, 0.0);
    JacobiSolver jac_solver(exacSol, params, argc, argv, rank, size);
    
    auto [sol, norm, n_it] = jac_solver.solve();
    
    if(rank == 0)
    {
        std::cout << "Norm: " << norm << std::endl;
        std::cout << "Number of iterations: " << n_it << std::endl;

        generateVTKFile("../files/output.vtk", sol, N, N, h, h);
    }

    MPI_Finalize();
    
    return 0;
}