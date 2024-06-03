#ifndef HH_INITIALIZEPROBLEM_HH
#define HH_INITIALIZEPROBLEM_HH

// this function initializes the problem defining an Eigen matrix of zeros except for the boundary

#include <Eigen/Dense>
#include "ParameterHandler.hpp"

void
initialize_problem(const param::ParameterHandler &params, Eigen::MatrixXd &M)
{
    const unsigned int n = params.coefficients.n; // number of intervals in the grid
    const double h = 1.0 / n; // step size
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
    std::array<double, 2> point;
    for(unsigned int j = 0; j < n; ++j)
    {
        point = {0.0, j * h};
        M(0, j) = params.functions.funBC_1(point); // bottom boundary
        
        point = {1.0, j * h};
        M(n - 1, j) = params.functions.funBC_3(point); // up boundary
    }
    for(unsigned int i = 0; i < n; ++i)
    {   
        point = {i * h, 0.0};
        M(i, 0) = params.functions.funBC_4(point); // left boundary

        point = {i * h, 1.0};
        M(i, n - 1) = params.functions.funBC_2(point); // right boundary
    }

}



#endif // HH_INITIALIZEPROBLEM_HH