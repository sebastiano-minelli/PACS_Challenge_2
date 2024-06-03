#include <iostream>
#include "ParameterHandler.hpp"
#include <Eigen/Dense>

// here we want to solve locally a subdomain of the problem

class Solve_local
{
    private:

        const std::size_t n_rows; // number of rows
        const std::size_t n_cols; // number of columns
        const param::ParameterHandler params; // parameters
        const Eigen::MatrixXd L_loc; // local matrix
        double norm_loc = std::numeric_limits<double>::infinity(); // local norm

    public:

        Solve_local(const Eigen::MatrixXd &L_loc_, const param::ParameterHandler &params_)
        : 
        n_rows(L_loc_.rows()), 
        n_cols(L_loc_.cols()), 
        params(params_), 
        L_loc(L_loc_) 
        {
            norm_loc = std::numeric_limits<double>::infinity();
        };
        
        // solves the laplace problem with Jacobi iteration locally
        std::pair<Eigen::MatrixXd, double> 
        solve()
        {
            const double h = 1.0 / params.coefficients.n; // step size

            Eigen::MatrixXd L_loc_new(n_rows, n_cols); // new local matrix

            for(std::size_t k = 0; k < params.coefficients.max_it; ++k) // if we haven't reached the max iterations
            {
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
                                        h * h * params.functions.fun({i, j})
                                        );
                        norm_loc += (L_loc_new(i, j) - L_loc(i, j)) * (L_loc_new(i, j) - L_loc(i, j));  
                        // notice that since we aren't changing the boundary its contribution to norm_loc is zero                  
                    }
                }
                norm_loc = std::sqrt(h * norm_loc);
                if(norm_loc < params.coefficients.tol_res)
                    break;
            }

            return std::pair<Eigen::MatrixXd, double>(L_loc_new, norm_loc);
        }







};