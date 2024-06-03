#include "ParameterHandler.hpp"
#include "writeVTK.hpp"
#include "Eigen/Dense"
#include <iostream>
#include "mpi_utils.hpp"
#include "partitioner.hpp"
#include "SafeMPI.hpp"
#include "InitializeProblem.hpp"
#include "LocalSolver.hpp"

int main()
{
    param::ParameterHandler params("data.txt");
    params.show_data();

    const int N = params.coefficients.n;

    const double h = 1.0/ N;

    Eigen::MatrixXd exacSol(N, N);
    initialize_problem(params, exacSol); // exacSol is passed by reference
    // std::cout << "Exact solution generated" << "\n" << exacSol << std::endl;
    LocalSolver loc_solver(exacSol, params);
    auto [sol, norm, n_it] = loc_solver.solve();
    std::cout << "Local norm: " << norm << std::endl;
    std::cout << "Number of iterations: " << n_it << std::endl;
    // std::cout << "Local solution: " << "\n" << sol << std::endl;

    generateVTKFile("../files/output.vtk", sol, N, N, h, h);


    return 0;
}