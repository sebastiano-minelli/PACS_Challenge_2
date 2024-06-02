#include "ParameterHandler.hpp"
#include "writeVTK.hpp"
#include <iostream>

int main()
{
    ParameterHandler param("data.txt");
    param.show_data();

    const int N = param.coefficients.n;

    const double h = 1.0/ (N + 1);

    std::vector<std::vector<double>> exacSol;

    double x, y;
    std::array<double, DIM> xvalue;
    double funvalue;

    for (int i = 0; i < N + 1; i++) 
    {
        std::vector<double> row;
        x = i * h;
        for (int j = 0; j < N + 1; j++) 
        {
            y = j * h;
            xvalue[0] = x;
            xvalue[1] = y;
            funvalue = param.functions.fun(xvalue);
            row.push_back(funvalue);
        }
        exacSol.push_back(row);
    }

    generateVTKFile("../files/output.vtk", exacSol, N, N, h, h);


    return 0;
}