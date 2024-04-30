#include <iostream>
#include <algorithm>
#include <iomanip>
#include "Matrix.hpp"
#include "chrono.hpp"

int main(int argc, char *argv[])
{
    using var_type = double; // value type
    constexpr algebra::StorageOrder ordering = algebra::StorageOrder::ROW_WISE; // storage method
    constexpr algebra::NormType norm_type1 = algebra::NormType::One; // norm type one
    constexpr algebra::NormType norm_type2 = algebra::NormType::Infinity; // norm type two
    constexpr algebra::NormType norm_type3 = algebra::NormType::Frobenius; // norm type three

    Timings::Chrono clock;

    algebra::Matrix<var_type, ordering> M;

    std::string filename = "lnsp_131.mtx"; // uncomment this line to test with 'lnsp_131.mtx' file
    // std::string filename = "matrix1.mtx"; // uncomment this line to test with 'matrix1.mtx' file

    M.parse_from_file(filename);
    std::cout << "Number of rows:            " << M.rows() << std::endl;
    std::cout << "Number of columns:         " << M.columns() << std::endl;
    std::cout << "Storing method:            " << (
                ordering == algebra::StorageOrder::ROW_WISE ? "'ROW-WISE'" : "'COLUMN-WISE'") << std::endl;
    std::cout << std::endl;
    

    //////// performing multiplication tests ///////////////////

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

    std::cout << "----- Multiplication performance test results -----" << std::endl;
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Uncompressed best time:    " << time_unc << " microseconds" << std::endl;
    std::cout << "Compressed best time:      " << time_comp << " microseconds" << std::endl;
    std::cout << std::endl;


    ///////// performing norm tests ////////////////////
    std::cout << "--------------- Norms test results ----------------" << std::endl;
    // One norm
    M.uncompress();
    var_type norm11 = M.norm<norm_type1>();
    std::cout << "Norm type:                 " << (
                norm_type1 == algebra::NormType::One ? "'ONE'" : 
                    norm_type1 == algebra::NormType::Infinity ? "'INFINITY'" : "'FROBENIUS'") << std::endl;
    std::cout << "Uncompressed state" << std::endl;
    std::cout << "Norm:                      " << norm11 << std::endl;

    M.compress();
    var_type norm12 = M.norm<norm_type1>();
    std::cout << "Compressed state" << std::endl;
    std::cout << "Norm:                      " << norm12 << std::endl;
    std::cout << std::endl;
    // Infinity norm
    M.uncompress();
    var_type norm21 = M.norm<norm_type2>();
    std::cout << "Norm type:                 " << (
                norm_type2 == algebra::NormType::One ? "'ONE'" : 
                    norm_type2 == algebra::NormType::Infinity ? "'INFINITY'" : "'FROBENIUS'") << std::endl;
    std::cout << "Uncompressed state" << std::endl;
    std::cout << "Norm:                      " << norm21 << std::endl;

    M.compress();
    var_type norm22 = M.norm<norm_type2>();
    std::cout << "Compressed state" << std::endl;
    std::cout << "Norm:                      " << norm22 << std::endl;
    std::cout << std::endl;
    // Frobenius norm
    M.uncompress();
    var_type norm31 = M.norm<norm_type3>();
    std::cout << "Norm type:                 " << (
                norm_type3 == algebra::NormType::One ? "'ONE'" : 
                    norm_type3 == algebra::NormType::Infinity ? "'INFINITY'" : "'FROBENIUS'") << std::endl;
    std::cout << "Uncompressed state" << std::endl;
    std::cout << "Norm:                      " << norm31 << std::endl;

    M.compress();
    var_type norm32 = M.norm<norm_type3>();
    std::cout << "Compressed state" << std::endl;
    std::cout << "Norm:                      " << norm32 << std::endl;
    std::cout << std::endl;
    
    return 0;
}