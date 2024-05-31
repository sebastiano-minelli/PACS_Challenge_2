#ifndef HH_PARAMETER_HANDLER_HH
#define HH_PARAMETER_HANDLER_HH

#include "Parameters.hpp"
#include<vector>
#include<array>

template<unsigned int DIM>
class ParameterHandler : public Parameters<DIM>
{
public:
    ParameterHandler(const std::string &filename) :
        Parameters<DIM>(filename)
    {};

    // function that handles automatically the analytic or numeric gradient
    const std::vector<double> compute_gradient(const std::array<double, DIM> point) const
    {
        std::vector<double> dfun_eval(DIM, 0.0);
        
        if(!this->coefficients.compute_num_grad)
        {
            for(size_t i = 0; i < DIM; ++i)
                dfun_eval[i] = this->function_param.dfun[i](point);
        }else
        {
            for(size_t i = 0; i < DIM; ++i)
            {
                std::array<double, DIM> x_plus = point;
                std::array<double, DIM> x_minus = point;
                x_plus[i] += this->coefficients.h;
                x_minus[i] -= this->coefficients.h;
                dfun_eval[i] = (this->function_param.fun(x_plus) - this->function_param.fun(x_minus)) /  (2 * this->coefficients.h);
            }
        }
        
        return dfun_eval;
    };

    void show_data() const
    {
        std::cout << "--------------- SELECTED OPTIONS ---------------\n" << std::endl;
        std::cout << "- Function:                     " << this->function_param.funString << "\n" << std::endl;
        if(!this->coefficients.compute_num_grad)
        {
            std::cout << "- Gradient of the function:     {1}: " << this->function_param.dfunString[0] << std::endl;
            for(size_t i = 1; i < DIM; ++i)
                std::cout<< "                                {" << i + 1 << "}: " << this->function_param.dfunString[i] << std::endl;
        }else
        {
            std::cout << "- Gradient of the function:     none" << std::endl;
        }
        std::cout << "\n" << std::endl;
        std::cout << "- Initial point:                {1}: " << this->function_param.x[0] << std::endl;
        for(size_t i = 1; i < DIM; ++i)
            std::cout<< "                                {" << i + 1 << "}: " << this->function_param.x[i] << std::endl;
        std::cout << "\n" << std::endl;
        std::cout << "- Maximum n. of iterations:     " << this->coefficients.max_it << "\n" << std::endl;
        std::cout << "- Residue tolerance:            " << this->coefficients.tol_res << "\n" << std::endl;
        std::cout << "- Argument tolerance (L2 norm): " << this->coefficients.tol_x << "\n" << std::endl;
        std::cout << "- Compute numeric gradient:     " << this->coefficients.compute_num_grad << "\n" << std::endl;
        std::cout << "- Finite difference increment:  " << this->coefficients.h << "\n" << std::endl;
        std::cout << "- Solver type:                  '" << this->min_method.solver_type << "'" << "\n" << std::endl;
        std::cout << "- Step coefficient method:      '" << this->step_coeff_method.coeff_solver << "'" << "\n" << std::endl;
        std::cout << "- Step coefficient alpha_0:     '" << this->step_coeff_method.alpha_zero << "'" << "\n" << std::endl;
        std::cout << "- Step coefficient mu:          '" << this->step_coeff_method.mu << "'" << "\n" << std::endl;
        std::cout << "- Step coefficient sigma:       '" << this->step_coeff_method.sigma << "'" << "\n" << std::endl;
        std::cout << "---------------------------------------------- \n" << std::endl;
    };
};



#endif