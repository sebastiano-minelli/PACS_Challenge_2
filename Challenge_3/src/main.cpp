#include "ParameterHandler.hpp"
#include "writeVTK.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include "mpi_utils.hpp"
#include "SafeMPI.hpp"
#include "JacobiSolver.hpp"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<double> L2_norms(5, 0.0); // vector for the norms (compared with the exact solution)
    std::vector<double> norms(5, 0.0); // vector for the norms
    std::vector<unsigned int> n_iterations(5, 0); // vector for the number of iterations
    std::vector<std::vector<double>> solutions; // vector for the matrices
    solutions.resize(5);

    for(int i = 0; i < 5; i++) 
    {
        std::string filename = "data" + std::to_string(i + 1) + ".txt";
        if(rank == 0)
            std::cout << "filename: "  << filename<< std::endl;

        param::ParameterHandler params(filename);
        /*
        if(rank == 0)
            params.show_data();
        */

        const int N = params.coefficients.n;

        const double h = 1.0/ N;

        std::vector<double> exact_solution = params.functions.fun_values; // exact solution
        JacobiSolver jac_solver(exact_solution, params, argc, argv, rank, size);
        
        std::tie(solutions[i], norms[i], n_iterations[i]) = jac_solver.solve();

        // calculate the L2 norm
        for(int x = 0; x < N; ++x)
        {
            for(int y = 0; y < N; ++y)
            {
                L2_norms[i] += (solutions[i][x * N + y] - exact_solution[x * N + y]) * 
                                (solutions[i][x * N + y] - exact_solution[x * N + y]);
            }
        }
        L2_norms[i] = std::sqrt(h * L2_norms[i]);
        
        if(rank == 0)
            generateVTKFile("../files/output" + std::to_string(i) + ".vtk", solutions[i], N, N, h, h);
    }

    if(rank == 0)
    {
        std::cout << "Test finished successfully\n" << std::endl;
        std::cout << std::setw(10) << "Matrix dimension" << std::setw(10) << "L2 norm" << std::endl;
        std::cout << std::endl;
        for(int i = 0; i < 5; i++)
        {
            std::cout << std::setw(10) << std::pow(2, i + 1)  << std::setw(16) << L2_norms[i] << std::endl;
        }
    }

    MPI_Finalize();
    
    return 0;
}