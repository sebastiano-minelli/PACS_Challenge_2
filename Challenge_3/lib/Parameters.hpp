#ifndef HH_PARAMETERS_HH
#define HH_PARAMETERS_HH

#include<iostream>
#include <string>
#include <vector>
#include "GetPot"
#include "muParserXInterface.hpp"
#include "Utilities.hpp"

using namespace MuParserInterface;

/*
BCs enumeration: 

     3
  ________
  |      |
4 |      | 2
  |      | 
  |______|
    1

*/

struct Function
{
    std::string funString; // force function string

    std::string funBC_1String; // BC of the bottom boundary string

    std::string funBC_2String; // BC of the right boundary string

    std::string funBC_3String; // BC of the up boundary string

    std::string funBC_4String; // BC of the left boundary string

    muParserXInterface<DIM> fun; // function

    muParserXInterface<DIM> funBC_1; // function

    muParserXInterface<DIM> funBC_2; // function

    muParserXInterface<DIM> funBC_3; // function

    muParserXInterface<DIM> funBC_4; // function

    std::array<double, DIM> x; // point of the domain
};

struct Coefficients
{
  unsigned int max_it; // maximum number of iterations
  
  double tol_res; // tolerance of the residual

  double tol_x; // tolerance of the argument

  double n; // number of interavals in the grid
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
    coefficients.h = datafile((section + "h").data(), 0.001);

    section = "Functions/";
    fun.funString = datafile((section + "fun").data(), " ");
    funBC_1.funString = datafile((section + "funBC_1").data(), " ");
    funBC_2.funString = datafile((section + "funBC_2").data(), " ");
    funBC_3.funString = datafile((section + "funBC_3").data(), " ");
    funBC_4.funString = datafile((section + "funBC_4").data(), " ");
    
    for(size_t i = 0; i < DIM; ++i)
    {
        std::string scalar_place = std::to_string(i + 1);
        function_param.x[i] = datafile((section + "x_" + scalar_place).data(), 0.0);
    }



    // Creating muParserX function and respective gradient
    MuParserInterface::muParserXInterface<DIM> dummy_fun(function_param.funString);
    function_param.fun = dummy_fun;
    
    if(!coefficients.compute_num_grad)
      for(size_t i = 0; i < DIM; ++i)
        function_param.dfun.emplace_back(function_param.dfunString[i]);
  }

  Function<DIM> fun;
  Function<DIM> funBC_1;
  Function<DIM> funBC_2;
  Function<DIM> funBC_3;
  Function<DIM> funBC_4;
  Coefficients coefficients;

#endif /* HH_PARAMETERS_HH */
