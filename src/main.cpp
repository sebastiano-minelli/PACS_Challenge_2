#include <iostream>
#include "DynamicMatrix.hpp"

int main()
{
    std::size_t n_rows = static_cast<std::size_t>(3); // numer of rows
    std::size_t n_cols = static_cast<std::size_t>(2); // number of columns
    std::map<std::array<std::size_t, 2>, double> elements; // elements

    constexpr algebra::StorageOrder order = algebra::ROW_WISE;

    elements[{0, 0}] = 1.0;
    elements[{0, 1}] = 2.0;
    elements[{0, 2}] = 3.0;
    elements[{1, 0}] = 4.0;
    elements[{1, 1}] = 5.0;
    elements[{1, 2}] = 6.0;
    elements[{2, 0}] = 7.0;
    elements[{2, 1}] = 8.0;
    elements[{2, 2}] = 9.0;

    algebra::DynamicMatrix<double, order> mat(n_rows, n_cols, std::move(elements));

    std::cout << "columns: " << mat.columns() << std::endl;
    std::cout << "rows: " << mat.rows() << std::endl;

    return 0;
}