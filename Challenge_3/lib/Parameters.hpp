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

struct Functions
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
    functions.funString = datafile((section + "fun").data(), "0.0 * x_1 * x_2");
    functions.funBC_1String = datafile((section + "funBC_1").data(), "0.0 * x_1 * x_2");
    functions.funBC_2String = datafile((section + "funBC_2").data(), "0.0 * x_1 * x_2");
    functions.funBC_3String = datafile((section + "funBC_3").data(), "0.0 * x_1 * x_2");
    functions.funBC_4String = datafile((section + "funBC_4").data(), "0.0 * x_1 * x_2");
  }

  Functions functions;
  Coefficients coefficients;
};

#endif /* HH_PARAMETERS_HH */
