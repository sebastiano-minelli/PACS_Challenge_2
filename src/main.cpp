#include <iostream>
#include <algorithm>
#include <iomanip>
#include "Matrix.hpp"
#include "chrono.hpp"

int main(int argc, char *argv[])
{
    using var_type = double; // value type
    constexpr algebra::StorageOrder ordering = algebra::StorageOrder::ROW_WISE; // storage method

    Timings::Chrono clock;

    algebra::Matrix<var_type, ordering> M;

    std::string filename = "lnsp_131.mtx";
    // std::string filename = "matrix1.mtx";

    M.parse_from_file(filename);
    std::cout << "Number of rows:       " << M.rows() << std::endl;
    std::cout << "Number of columns:    " << M.columns() << std::endl;
    

    // perform multiplication tests

    std::vector<var_type> v(M.columns(), static_cast<var_type>(1.0)); // vector of ones
    std::vector<var_type> res1(M.rows(), static_cast<var_type>(0.0)); // uncompressed result
    std::vector<var_type> res2(M.rows(), static_cast<var_type>(0.0)); // compressed result

    std::size_t n = 10000; // number of tests

    // uncompress state test
    double time_unc = std::numeric_limits<double>::max();
    M.uncompress();
    for(std::size_t i = 0; i < n; ++i)
    {
        clock.start();
        res2 = std::move(M * v);
        clock.stop();
        auto walltime = clock.wallTime();
        time_unc = std::min(walltime, time_unc);
    }

    // compress state test
    double time_comp = std::numeric_limits<double>::max();
    M.compress();
    for(std::size_t i = 0; i < n; ++i)
    {
        clock.start();
        res1 = std::move(M * v);
        clock.stop();
        auto walltime = clock.wallTime();
        time_comp = std::min(walltime, time_comp);
    }

    std::cout << "Multiplication performance test results:" << std::endl;
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Uncompressed best time:    " << time_unc << " microseconds" << std::endl;
    std::cout << "Compressed best time:      " << time_comp << " microseconds" << std::endl;
    
    
    return 0;
}