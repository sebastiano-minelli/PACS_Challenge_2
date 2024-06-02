#include "ParameterHandler.hpp"
#include <iostream>

int main()
{
    ParameterHandler param("data.txt");
    param.show_data();

    return 0;
}