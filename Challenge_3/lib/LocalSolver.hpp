#ifndef HH_LOCAL_SOLVER_HH
#define HH_LOCAL_SOLVER_HH

#include <iostream>
#include <vector>
#include "ParameterHandler.hpp"
#include <omp.h>

// here we want to solve locally a subdomain of the problem

class LocalSolver
{
    private:

        std::size_t n_rows; // number of rows
        std::size_t n_cols; // number of columns
        const int global_row_index; // idex of the first row of the global matrix owned by the process 
        const param::ParameterHandler params; // parameters
        const std::vector<double> L; // matrix
        double norm_loc = std::numeric_limits<double>::infinity(); // local norm

    public:

        LocalSolver(const std::vector<double> &L_, const param::ParameterHandler &params_, const int global_row_index_)
        :
        global_row_index(global_row_index_),
        params(params_), 
        L(L_),
        norm_loc(std::numeric_limits<double>::infinity())
        {
            // compute the number of rows and columns
            n_rows = L.size() / params.coefficients.n;
            n_cols = params.coefficients.n;
        };
        
        // solves the laplace problem with Jacobi iteration locally 
        std::tuple<std::vector<double>, double> 
        solve()
        {
            const double h = 1.0 / params.coefficients.n; // step size
            std::vector<double> L_loc = L; // local matrix
            std::vector<double> L_loc_new = L; // new local matrices

            
            norm_loc = 0.0; // reset the norm

            // making the loop parallel
            #pragma omp parallel for shared(L_loc_new) shared(L_loc) reduction(+:norm_loc)
            // compute Jacobi iteration
            for(std::size_t i = 1; i < n_rows - 1; ++i)
            {
                // Notice that muParserXInterface isn't thread safe, we have to define a parser in every thread
                MuParserInterface::muParserXInterface<2> parser(params.functions.fun);
                
                for(std::size_t j = 1; j < n_cols - 1; ++j)
                {
                    
                    std::array<double, 2> vars = {j * h, (global_row_index + i) * h};
                    L_loc_new[i * n_cols + j] = 0.25 * (
                                    L_loc[(i - 1) * n_cols + j] + 
                                    L_loc[(i + 1) * n_cols + j] + 
                                    L_loc[i * n_cols + j - 1] + 
                                    L_loc[i * n_cols + j + 1] +
                                    h * h * parser(vars) // notice that i and j are inverted
                                    );
                    norm_loc += (L_loc_new[i * n_cols + j] - L_loc[i * n_cols + j]) * (L_loc_new[i * n_cols + j] - L_loc[i * n_cols + j]);  
                    // notice that since we aren't changing the boundary its contribution to norm_loc is zero                  
                }
            }
            norm_loc = std::sqrt(h * norm_loc);
            // update
            L_loc = L_loc_new;
            

            return std::make_tuple(L_loc_new, norm_loc);
        }
};

#endif // HH_LOCAL_SOLVER_HH