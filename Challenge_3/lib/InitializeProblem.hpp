#ifndef HH_INITIALIZEPROBLEM_HH
#define HH_INITIALIZEPROBLEM_HH

// this function initializes the problem defining an Eigen matrix of zeros except for the boundary

#include <Eigen/Dense>
#include "ParameterHandler.hpp"

void
initialize_problem(const param::ParameterHandler &params, Eigen::MatrixXd &M)
{
    unsigned int n = params.coefficients.n; // number of intervals in the grid
    if(n != M.rows() || n != M.cols())
    {
        std::cerr << "Error: the matrix M must have the same number of rows as the number of intervals in the grid" << std::endl;
        std::exit(1);
    }
    if(n < 3)
    {
        std::cerr << "Error: the number of intervals in the grid must be at least 3" << std::endl;
        std::exit(1);
    }

    M.setZero(); 
    for(unsigned int j = 0; j < n; ++j)
    {
        M(0, j) = params.functions.funBC_1({0, j}); // bottom boundary
        M(n - 1, j) = params.functions.funBC_3({n - 1, j}); // up boundary
    }
    for(unsigned int i = 0; i < n; ++i)
    {
        M(i, 0) = params.functions.funBC_4({i, 0}); // left boundary
        M(i, n - 1) = params.functions.funBC_2({i, n - 1}); // right boundary
    }

}



#endif // HH_INITIALIZEPROBLEM_HH