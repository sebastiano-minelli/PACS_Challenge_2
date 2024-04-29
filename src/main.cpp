#include <iostream>
#include "Matrix.hpp"

int main(int argc, char *argv[])
{
    
    algebra::Matrix<double, algebra::StorageOrder::ROW_WISE> m_row;

    // std::string filename = "lnsp_131.mtx";
    std::string filename = "matrix1.mtx";

    m_row.parse_from_file(filename);

    std::cout << "Matrix read from file: " << std::endl;
    std::cout << m_row(0, 1) << std::endl;
    std::cout << m_row(0, 2) << std::endl;
    std::cout << m_row(0, 3) << std::endl;
    std::cout << m_row(1, 0) << std::endl;
    std::cout << m_row(1, 1) << std::endl;
    std::cout << m_row(2, 2) << std::endl;
    std::cout << "Compressed state: " << m_row.is_compressed() << std::endl;
    m_row.compress();
    std::cout << "Compressed state: " << m_row.is_compressed() << std::endl;
    std::cout << "Compressed matrix read from file: " << std::endl;
    std::cout << m_row(0, 1) << std::endl;
    std::cout << m_row(0, 2) << std::endl;
    std::cout << m_row(0, 3) << std::endl;
    std::cout << m_row(1, 0) << std::endl;
    std::cout << m_row(1, 1) << std::endl;
    std::cout << m_row(2, 2) << std::endl;
    


    m_row.compress();
    std::cout << "Compressed state: " << m_row.is_compressed() << std::endl;
    std::cout << "Uncompressed matrix read from file: " << std::endl;
    std::cout << m_row(0, 1) << std::endl;
    std::cout << m_row(0, 2) << std::endl;
    std::cout << m_row(0, 3) << std::endl;
    std::cout << m_row(1, 0) << std::endl;
    std::cout << m_row(1, 1) << std::endl;
    std::cout << m_row(2, 2) << std::endl;
    
    m_row.compress();
    std::cout << "Matrix-vector product: " << std::endl;
    std::vector<double> x = {2.0, 0.0, 1.0, 0.0};
    std::vector<double> y = m_row * x;
    std::cout << y[0] << std::endl;
    std::cout << y[1] << std::endl;
    std::cout << y[2] << std::endl;
    std::cout << y[3] << std::endl;

    std::cout << "Number of rows: " << m_row.rows() << std::endl;
    std::cout << "Number of columns: " << m_row.columns() << std::endl;

    // resizing test
    m_row.uncompress();
    m_row.resize(10, 10);
    std::cout << "Number of rows: " << m_row.rows() << std::endl;
    std::cout << "Number of columns: " << m_row.columns() << std::endl;
    std::cout << "Matrix read: " << std::endl;
    std::cout << m_row(0, 1) << std::endl;
    std::cout << m_row(0, 2) << std::endl;
    std::cout << m_row(0, 3) << std::endl;
    std::cout << m_row(1, 0) << std::endl;
    std::cout << m_row(1, 1) << std::endl;
    std::cout << m_row(2, 2) << std::endl;
    std::cout << m_row(9, 9) << std::endl;

    



    return 0;
}