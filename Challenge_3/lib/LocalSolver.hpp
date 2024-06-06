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

        const int n_rows; // number of rows
        const int n_cols; // number of columns
        const int global_row_index; // idex of the first row of the global matrix owned by the process 
        const param::ParameterHandler params; // parameters
        double norm_loc = std::numeric_limits<double>::infinity(); // local norm

    public:

        LocalSolver(const int n_rows_, const int n_cols_, const param::ParameterHandler &params_, const int global_row_index_)
        :
        n_rows(n_rows_),
        n_cols(n_cols_),
        global_row_index(global_row_index_),
        params(params_),
        norm_loc(std::numeric_limits<double>::infinity())
        {};
        
        // solves the laplace problem with Jacobi iteration locally 
        double solve(std::vector<double> &L_loc)
        {
            const double h = 1.0 / params.coefficients.n; // step size
            std::vector<double> L_loc_new = L_loc; // new local matrices

            norm_loc = 0.0; // reset the norm

            // making the loop parallel
            #pragma omp parallel for shared(L_loc_new) shared(L_loc) reduction(+:norm_loc)
            // compute Jacobi iteration
            for(int i = 1; i < n_rows - 1; ++i)
            {                
                for(int j = 1; j < n_cols - 1; ++j)
                {
                    L_loc_new[i * n_cols + j] = 0.25 * (
                                    L_loc[(i - 1) * n_cols + j] + 
                                    L_loc[(i + 1) * n_cols + j] + 
                                    L_loc[i * n_cols + j - 1] + 
                                    L_loc[i * n_cols + j + 1] +
                                    h * h * params.functions.fun_values[(global_row_index + i) * n_cols + j + 1]
                                    );
                    norm_loc += (L_loc_new[i * n_cols + j] - L_loc[i * n_cols + j]) * (L_loc_new[i * n_cols + j] - L_loc[i * n_cols + j]);  
                    // notice that since we aren't changing the boundary its contribution to norm_loc is zero                  
                }
            }
            norm_loc = std::sqrt(h * norm_loc);

            // update
            L_loc = L_loc_new;

            return norm_loc;
        }
};

#endif // HH_LOCAL_SOLVER_HH