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
        const int global_row_index; // index of the first row of the global matrix owned by the process 
        const param::ParameterHandler params; // parameters
        double local_norm = std::numeric_limits<double>::infinity(); // local norm

    public:

        LocalSolver(const int n_rows_, const int n_cols_, const param::ParameterHandler &params_, const int global_row_index_)
        :
        n_rows(n_rows_),
        n_cols(n_cols_),
        global_row_index(global_row_index_),
        params(params_),
        local_norm(std::numeric_limits<double>::infinity())
        {};
        
        // solves the laplace problem with Jacobi iteration locally (returns the new local matrix and the local norm)
        std::tuple<std::vector<double>, double> 
        solve(const std::vector<double> &L_loc)
        {
            const double h = 1.0 / params.coefficients.n; // step size
            std::vector<double> L_loc_new = L_loc; // new local matrices

            local_norm = 0.0; // reset the norm

            // making the loop parallel
            #pragma omp parallel for shared(L_loc_new) shared(L_loc) reduction(+:local_norm)
            for(int i = 1; i < n_rows - 1; ++i)
            {                
                for(int j = 1; j < n_cols - 1; ++j)
                {
                    L_loc_new[i * n_cols + j] = 0.25 * 
                                                (
                                                L_loc[(i - 1) * n_cols + j] + 
                                                L_loc[(i + 1) * n_cols + j] + 
                                                L_loc[i * n_cols + j - 1] + 
                                                L_loc[i * n_cols + j + 1] +
                                                h * h * params.functions.fun_values[(global_row_index + i) * n_cols + j]
                                                );
                    local_norm += (L_loc_new[i * n_cols + j] - L_loc[i * n_cols + j]) * (L_loc_new[i * n_cols + j] - L_loc[i * n_cols + j]) / 
                                  (L_loc[i * n_cols + j] * L_loc[i * n_cols + j]);  
                    // notice that since we aren't changing the boundary its contribution to norm_loc is zero                  
                }
            }
            
            local_norm = std::sqrt(h * local_norm);

            return std::make_tuple(L_loc_new, local_norm);
        }
};

#endif // HH_LOCAL_SOLVER_HH