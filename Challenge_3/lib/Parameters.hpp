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
    Domain enumeration
        ----3----
        |       |
        4       2
        |       |
        ----1----
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

    std::vector<double> funBC_1_values; // values of the BC1

    std::vector<double> funBC_2_values; // values of the BC2

    std::vector<double> funBC_3_values; // values of the BC3

    std::vector<double> funBC_4_values; // values of the BC4

    std::vector<double> fun_values; // values of the function
};

struct Coefficients
{
  unsigned int max_it; // maximum number of iterations
  
  double tol_res; // tolerance of the residual

  unsigned int n; // number of intervals in the grid
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
    coefficients.n = datafile((section + "n").data(), 10);

    section = "Functions/";
    functions.funString = datafile((section + "fun").data(), "0.0 * x[1] * x[2]");
    functions.funBC_1String = datafile((section + "funBC_1").data(), "0.0 * x[1] * x[2]");
    functions.funBC_2String = datafile((section + "funBC_2").data(), "0.0 * x[1] * x[2]");
    functions.funBC_3String = datafile((section + "funBC_3").data(), "0.0 * x[1] * x[2]");
    functions.funBC_4String = datafile((section + "funBC_4").data(), "0.0 * x[1] * x[2]");

    // Creating muParserX function and respective gradient
    MuParserInterface::muParserXInterface<DIM> dummy_fun(functions.funString);
    functions.fun = dummy_fun;

    // Evaluating the function
    functions.fun_values.resize(coefficients.n * coefficients.n);
    unsigned int n = coefficients.n;
    double h = 1.0 / n;
    std::array<double, 2> vars;
    for (unsigned int i = 0; i < n; ++i)
    {
      for (unsigned int j = 0; j < n; ++j)
      {
        vars = {j * h, i * h};
        functions.fun_values[i * n + j] = functions.fun(vars);
      }
    }

    MuParserInterface::muParserXInterface<DIM> dummy_fun_BC1(functions.funBC_1String);
    functions.funBC_1 = dummy_fun_BC1;

    // Evaluating the BC1
    functions.funBC_1_values.resize(coefficients.n);
    for (unsigned int i = 0; i < n; ++i)
    {
      vars = {i * h, 0.0};
      functions.funBC_1_values[i] = functions.funBC_1(vars);
    }

    MuParserInterface::muParserXInterface<DIM> dummy_fun_BC2(functions.funBC_2String);
    functions.funBC_2 = dummy_fun_BC2;

    // Evaluating the BC2
    functions.funBC_2_values.resize(coefficients.n);
    for (unsigned int i = 0; i < n; ++i)
    {
      vars = {1.0, i * h};
      functions.funBC_2_values[i] = functions.funBC_2(vars);
    }

    MuParserInterface::muParserXInterface<DIM> dummy_fun_BC3(functions.funBC_3String);
    functions.funBC_3 = dummy_fun_BC3;

    // Evaluating the BC3
    functions.funBC_3_values.resize(coefficients.n);
    for (unsigned int i = 0; i < n; ++i)
    {
      vars = {i * h, 1.0};
      functions.funBC_3_values[i] = functions.funBC_2(vars);
    }

    MuParserInterface::muParserXInterface<DIM> dummy_fun_BC4(functions.funBC_4String);
    functions.funBC_4 = dummy_fun_BC4;

    // Evaluating the BC4
    functions.funBC_4_values.resize(coefficients.n);
    for (unsigned int i = 0; i < n; ++i)
    {
      vars = {0.0, i * h};
      functions.funBC_4_values[i] = functions.funBC_2(vars);
    }
  }

  Functions functions;
  Coefficients coefficients;
};

#endif /* HH_PARAMETERS_HH */
