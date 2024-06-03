#ifndef HH_PARAMETER_HANDLER_HH
#define HH_PARAMETER_HANDLER_HH

#include "Parameters.hpp"
#include<vector>
#include<array>

namespace param
{
class ParameterHandler : public Parameters
{
public:
    ParameterHandler(const std::string &filename) : Parameters(filename)
    {};

    void show_data() const
    {
        std::cout << "--------------- SELECTED OPTIONS ---------------\n" << std::endl;
        std::cout << "- Function:                     " << this->functions.funString << "\n" << std::endl;
        std::cout << "Domain boundaries enumeration" << std::endl;
        std::cout << "      ----3----" << std::endl;
        std::cout << "      |       |" << std::endl;
        std::cout << "      4       2" << std::endl;
        std::cout << "      |       |" << std::endl;
        std::cout << "      ----1----" << "\n" << std::endl;
        std::cout << "- BC on 1:                      " << this->functions.funBC_1String << "\n" << std::endl;
        std::cout << "- BC on 2:                      " << this->functions.funBC_2String << "\n" << std::endl;
        std::cout << "- BC on 3:                      " << this->functions.funBC_3String << "\n" << std::endl;
        std::cout << "- BC on 4:                      " << this->functions.funBC_4String << "\n" << std::endl;

        std::cout << "\n" << std::endl;
        std::cout << "\n" << std::endl;
        std::cout << "- Maximum n. of iterations:     " << this->coefficients.max_it << "\n" << std::endl;
        std::cout << "- Residue tolerance:            " << this->coefficients.tol_res << "\n" << std::endl;
        std::cout << "- Argument tolerance (L2 norm): " << this->coefficients.tol_x << "\n" << std::endl;
        std::cout << "- Number of grid intervals:     " << this->coefficients.n << "\n" << std::endl;
        std::cout << "---------------------------------------------- \n" << std::endl;
    };
};

} // end namespace param



#endif