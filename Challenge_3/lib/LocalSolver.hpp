#ifndef HH_LOCAL_SOLVER_HH
#define HH_LOCAL_SOLVER_HH

#include <iostream>
#include "ParameterHandler.hpp"
#include <Eigen/Dense>

// here we want to solve locally a subdomain of the problem

class LocalSolver
{
    private:

        const std::size_t n_rows; // number of rows
        const std::size_t n_cols; // number of columns
        const param::ParameterHandler params; // parameters
        const Eigen::MatrixXd L; // matrix
        double norm_loc = std::numeric_limits<double>::infinity(); // local norm

    public:

        LocalSolver(const Eigen::MatrixXd &L_, const param::ParameterHandler &params_)
        : 
        n_rows(L_.rows()), 
        n_cols(L_.cols()), 
        params(params_), 
        L(L_) 

        {
            norm_loc = std::numeric_limits<double>::infinity();
        };
        
        // solves the laplace problem with Jacobi iteration locally 
        std::tuple<Eigen::MatrixXd, double, unsigned int> 
        solve()
        {
            const double h = 1.0 / params.coefficients.n; // step size

            unsigned int n_it = 0; // number of iterations

            Eigen::MatrixXd L_loc = L; // local matrix
            Eigen::MatrixXd L_loc_new = L; // new local matrices

            for(std::size_t k = 0; k < params.coefficients.max_it && norm_loc > params.coefficients.tol_res; ++k) // if we haven't reached the max iterations
            {
                norm_loc = 0.0; // reset the norm
                // compute Jacobi iteration
                for(std::size_t i = 1; i < n_rows - 1; ++i)
                {
                    for(std::size_t j = 1; j < n_cols - 1; ++j)
                    {
                        L_loc_new(i, j) = 0.25 * (
                                        L_loc(i - 1, j) + 
                                        L_loc(i + 1, j) + 
                                        L_loc(i, j - 1) + 
                                        L_loc(i, j + 1) +
                                        h * h * params.functions.fun({j * h, i * h}) // notice that i and j are inverted
                                        );
                        norm_loc += (L_loc_new(i, j) - L_loc(i, j)) * (L_loc_new(i, j) - L_loc(i, j));  
                        // notice that since we aren't changing the boundary its contribution to norm_loc is zero                  
                    }
                }
                norm_loc = std::sqrt(h * norm_loc);
                ++n_it;
                // update
                L_loc = L_loc_new;
            }

            return std::make_tuple(L_loc_new, norm_loc, n_it);
        }
};

#endif // HH_LOCAL_SOLVER_HH