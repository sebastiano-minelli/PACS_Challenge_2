#include <iostream>
#include "Matrix.hpp"

int main()
{
    algebra::Matrix<double, algebra::StorageOrder::ROW_WISE> m_row;

    // std::string filename = "lnsp_131.mtx";
    std::string filename = "matrix1.mtx";

    m_row.parse_from_file(filename);

    std::cout << "Matrix read from file: " << std::endl;
    std::cout << m_row(0, 0) << std::endl;
    std::cout << m_row(0, 1) << std::endl;
    std::cout << m_row(1, 0) << std::endl;
    std::cout << m_row(1, 1) << std::endl;
    std::cout << "Compressed state: " << m_row.is_compressed() << std::endl;
    m_row.compress();
    std::cout << "Compressed state: " << m_row.is_compressed() << std::endl;

    std::vector<double> x = {2.0, 0.0, 1.0};
    std::vector<double> y = m_row * x;
    std::cout << "Matrix-vector product: " << std::endl;
    std::cout << y[0] << std::endl;
    std::cout << y[1] << std::endl;
    std::cout << y[2] << std::endl;

    return 0;
}