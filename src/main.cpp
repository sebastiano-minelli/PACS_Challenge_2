#include <iostream>
#include "DynamicMatrix.hpp"
#include <iostream>
#include<map>
#include<array>

int main()
{
    std::size_t n_rows; // numer of rows
    std::size_t n_cols; // number of columns
    std::map<std::array<std::size_t, 2>, double> elements; // elements

    elements[{0, 0}] = 1.0;
    elements[{0, 1}] = 2.0;
    elements[{0, 2}] = 3.0;
    elements[{1, 0}] = 4.0;
    elements[{1, 1}] = 5.0;
    elements[{1, 2}] = 6.0;
    elements[{2, 0}] = 7.0;
    elements[{2, 1}] = 8.0;
    elements[{2, 2}] = 9.0;

    return 0;
}