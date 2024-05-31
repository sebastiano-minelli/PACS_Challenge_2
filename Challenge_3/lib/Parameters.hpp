#ifndef HH_PARAMETERS_HH
#define HH_PARAMETERS_HH

#include<iostream>
#include <string>
#include <vector>
#include "GetPot"
#include "muParserXInterface.hpp"

using namespace MuParserInterface;

template<unsigned int DIM>
struct Function
{
    std::string funString; // function string

    std::vector<std::string> dfunString; // gradient of the function string

    muParserXInterface<DIM> fun; // function

    std::vector<muParserXInterface<DIM>> dfun; // function gradient

    std::array<double, DIM> x; // point of the domain
};

struct Coefficients
{
  unsigned int max_it; // maximum number of iterations
  
  double tol_res; // tolerance of the residual

  double tol_x; // tolerance of the argument

  bool compute_num_grad; // true to compute a numerical gradient, false otherwise

  double h; // spacing for central difference scheme
};

struct StepCoeffMethod
{ // method to use to solve the step coefficient parameter
  std::string coeff_solver; // 'Exponential', 'InverseDecay', 'Armijo' 
  double alpha_zero;
  double mu;
  double sigma;
};

struct MinMethod
{ // method to use to solve the minimization problem
  std::string solver_type; // 'Gradient' 
};

template<unsigned int DIM>
class Parameters
{
public:
  Parameters(const std::string &filename)
  {
    // GetPot reads from file
    GetPot datafile(filename.c_str());

    std::string section = "Parameters/";

    coefficients.max_it = datafile((section + "max_it").data(), 500);
    coefficients.tol_res = datafile((section + "tol_res").data(), 1.0e-5);
    coefficients.tol_x = datafile((section + "tol_x").data(), 1.0e-5);
    coefficients.compute_num_grad = datafile((section + "compute_num_grad").data(), true);
    coefficients.h = datafile((section + "h").data(), 0.001);

    section = "Functions/";
    function_param.funString = datafile((section + "fun").data(), " ");
    
    if(!coefficients.compute_num_grad)
    {  
      function_param.dfunString.resize(DIM);
      for(size_t i = 0; i < DIM; ++i)
      {
          // evaluate the right term of the gradient vector
          std::string scalar_place = std::to_string(i + 1);
          function_param.dfunString[i] = datafile((section + "grad_fun_" + scalar_place).data(), " ");
      }
    }
    
    for(size_t i = 0; i < DIM; ++i)
    {
        std::string scalar_place = std::to_string(i + 1);
        function_param.x[i] = datafile((section + "x_" + scalar_place).data(), 0.0);
    }
    

    section = "Solver_Type/";
    min_method.solver_type = datafile((section + "solver_type").data(), "Gradient"); // Solver choice
 
    section = "Compute_Coefficient_Method/";
    step_coeff_method.coeff_solver = datafile((section + "coeff_solver").data(), "Exponential"); // Step coefficient choice
    step_coeff_method.alpha_zero = datafile((section + "Method_Options/" + "alpha_zero").data(), 0.5); //alpha zero
    step_coeff_method.mu = datafile((section + "Method_Options/" + "mu").data(), 0.2); //mu
    step_coeff_method.sigma = datafile((section + "Method_Options/" + "sigma").data(), 0.3); //sigma



    // Creating muParserX function and respective gradient
    MuParserInterface::muParserXInterface<DIM> dummy_fun(function_param.funString);
    function_param.fun = dummy_fun;
    
    if(!coefficients.compute_num_grad)
      for(size_t i = 0; i < DIM; ++i)
        function_param.dfun.emplace_back(function_param.dfunString[i]);
  }

  Function<DIM> function_param;
  Coefficients coefficients;
  StepCoeffMethod step_coeff_method;
  MinMethod min_method;
};

#endif /* HH_PARAMETERS_HH */
