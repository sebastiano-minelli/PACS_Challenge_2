#ifndef HH_INITIALIZEPROBLEM_HH
#define HH_INITIALIZEPROBLEM_HH

// this function initializes the problem defining an Eigen matrix of zeros except for the boundary

#include <vector>
#include "ParameterHandler.hpp"

void
initialize_problem(const param::ParameterHandler &params, std::vector<double> &M)
{
    const unsigned int n = params.coefficients.n; // number of intervals in the grid  

    std::fill(M.begin(), M.end(), 0.0); // fill the matrix with zeros

    #pragma omp parallel for shared(M)
    for(unsigned int i = 0; i < n; ++i)
    {
        M[n*n - n + i] = params.functions.funBC_1_values[i]; // bottom boundary // n*n - n + j = n(n-1) + j
        
        M[i] = params.functions.funBC_3_values[i]; // up boundary // 0*n + i = i

        M[i * n + n - 1] = params.functions.funBC_4_values[i]; // left boundary

        M[i * n] = params.functions.funBC_2_values[i]; // right boundary
    }
}

#endif // HH_INITIALIZEPROBLEM_HH