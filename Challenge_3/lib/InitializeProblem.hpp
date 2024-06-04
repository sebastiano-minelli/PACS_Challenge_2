#ifndef HH_INITIALIZEPROBLEM_HH
#define HH_INITIALIZEPROBLEM_HH

// this function initializes the problem defining an Eigen matrix of zeros except for the boundary

#include <Eigen/Dense>
#include "ParameterHandler.hpp"

void
initialize_problem(const param::ParameterHandler &params, Eigen::MatrixXd &M, int row_start, int row_end, int col_start, int col_end)
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
    #pragma omp parallel for shared(M)
    for(unsigned int j = col_start; j < col_end; ++j)
    {
        // notice that muParserXInterface isn't thread safe, we have to define a parser in every thread
        MuParserInterface::muParserXInterface<2> BC1(params.functions.funBC_1);
        MuParserInterface::muParserXInterface<2> BC3(params.functions.funBC_3);

        point = {j * h, 0.0};
        M(n - 1, j) = BC1(point); // bottom boundary
        
        point = {j * h, 1.0};
        M(0, j) = BC3(point); // up boundary
    }
    for(unsigned int i = row_start; i < row_end; ++i)
    {   
        // notice that muParserXInterface isn't thread safe, we have to define a parser in every thread
        MuParserInterface::muParserXInterface<2> BC2(params.functions.funBC_2);
        MuParserInterface::muParserXInterface<2> BC4(params.functions.funBC_4);

        point = {0.0, i * h};
        M(i, 0) = BC4(point); // left boundary

        point = {1.0, i * h};
        M(i, n - 1) = BC2(point); // right boundary
    }

}



#endif // HH_INITIALIZEPROBLEM_HH