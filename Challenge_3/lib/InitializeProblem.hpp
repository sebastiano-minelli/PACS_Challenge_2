#ifndef HH_INITIALIZEPROBLEM_HH
#define HH_INITIALIZEPROBLEM_HH

// this function initializes the problem defining an Eigen matrix of zeros except for the boundary

#include <vector>
#include "ParameterHandler.hpp"

void
initialize_problem(const param::ParameterHandler &params, std::vector<double> &M)
{
    const unsigned int n = params.coefficients.n; // number of intervals in the grid
    const double h = 1.0 / n; // step size    

    std::fill(M.begin(), M.end(), 0.0); // fill the matrix with zeros

    std::array<double, 2> point;
    #pragma omp parallel for shared(M)
    for(unsigned int i = 0; i < n; ++i)
    {
        // notice that muParserXInterface isn't thread safe, we have to define a parser in every thread
        MuParserInterface::muParserXInterface<2> BC1(params.functions.funBC_1);
        MuParserInterface::muParserXInterface<2> BC3(params.functions.funBC_3);
        MuParserInterface::muParserXInterface<2> BC2(params.functions.funBC_2);
        MuParserInterface::muParserXInterface<2> BC4(params.functions.funBC_4);

        point = {i * h, 0.0};
        M[n*n - n + i] = BC1(point); // bottom boundary // n*n - n + j = n(n-1) + j
        
        point = {i * h, 1.0};
        M[i] = BC3(point); // up boundary // 0*n + i = i

        point = {0.0, i * h};
        M[i * n + n - 1] = BC4(point); // left boundary

        point = {1.0, i * h};
        M[i * n] = BC2(point); // right boundary
    }

}



#endif // HH_INITIALIZEPROBLEM_HH