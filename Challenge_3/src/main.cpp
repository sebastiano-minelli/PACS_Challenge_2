#include "ParameterHandler.hpp"
#include "writeVTK.hpp"
#include "Eigen/Dense"
#include <iostream>
#include "mpi_utils.hpp"
#include "partitioner.hpp"
#include "SafeMPI.hpp"

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

    Eigen::MatrixXd L(N, N);

    apsc::MatrixPartitioner mpartitioner(3, 3, 2);
    std::cout << "First row: " << mpartitioner.first_row(5) << std::endl;


    return 0;
}