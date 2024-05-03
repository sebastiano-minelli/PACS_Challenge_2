#include <iostream>
#include <algorithm>
#include <iomanip>
#include "Matrix.hpp"
#include "chrono.hpp"

int main(int argc, char *argv[])
{
    ////////////////////// TEMPLATE VARIABLES AND FILES SELECTION ////////////////////////////////

    using var_type = double; // value type

    // parsing test
    std::string filename_parse = "../files/lnsp_131.mtx"; // matrix for the parsing test
    constexpr algebra::StorageOrder ordering_parse = algebra::StorageOrder::ROW_WISE; // storage method

    // Multiplication performance test matrix by vector
    std::string filename_performance = "../files/lnsp_131.mtx"; // matrix for the multiplication performance test
    constexpr algebra::StorageOrder ordering_performance = algebra::StorageOrder::ROW_WISE; // storage method

    // multiplication test matrix by vector with operator*(Matrix M, Matrix v)
    std::string filename_mat = "../files/matrix1.mtx"; // matrix for the multiplication test (Matrix -Matrix type)
    std::string filename_vector = "../files/vector_sparse.mtx"; // sparse vector for the matrix multiplication
    constexpr algebra::StorageOrder ordering_mat = algebra::StorageOrder::ROW_WISE; // storage method for the matrix
    constexpr algebra::StorageOrder ordering_vector = algebra::StorageOrder::ROW_WISE; // storage method for the vector

    // Norm test
    std::string filename_norms = "../files/lnsp_131.mtx"; // matrix for the norm test
    constexpr algebra::StorageOrder ordering_norms = algebra::StorageOrder::ROW_WISE; // storage method for the matrix
    constexpr algebra::NormType norm_type1 = algebra::NormType::One; // norm type one
    constexpr algebra::NormType norm_type2 = algebra::NormType::Infinity; // norm type two
    constexpr algebra::NormType norm_type3 = algebra::NormType::Frobenius; // norm type three    

    ///////////////////////////////////////////////////////////////////////////////////////////////


    ////////// performing parsing test ///////////////////
    algebra::Matrix<var_type, ordering_parse> M_parse;
    M_parse.parse_from_file(filename_parse);

    // output
    std::cout << "---------- Parsing test results ----------" << std::endl;
    std::cout << "Matrix read from file:     '" << filename_parse << "'" << std::endl;
    std::cout << "Number of rows:            " << M_parse.rows() << std::endl;
    std::cout << "Number of columns:         " << M_parse.columns() << std::endl;
    std::cout << "Storing method:            " << (
                ordering_parse == algebra::StorageOrder::ROW_WISE ? "'ROW-WISE'" : "'COLUMN-WISE'") << std::endl;
    std::cout << std::endl;
    

    //////// performing multiplication tests (matrix by vector) ///////////////////
    Timings::Chrono clock; // setting clock
    algebra::Matrix<var_type, ordering_performance> M_performance; // creating matrix
    M_performance.parse_from_file(filename_performance); // parsing matrix

    std::vector<var_type> v(M_performance.columns(), static_cast<var_type>(1.0)); // vector of ones
    std::vector<var_type> res1(M_performance.rows(), static_cast<var_type>(0.0)); // uncompressed result
    std::vector<var_type> res2(M_performance.rows(), static_cast<var_type>(0.0)); // compressed result

    std::size_t n = 10000; // number of tests

    // uncompress state test
    double time_unc = std::numeric_limits<double>::max();
    M_performance.uncompress();
    for(std::size_t i = 0; i < n; ++i)
    {
        clock.start();
        res2 = std::move(M_performance * v);
        clock.stop();
        auto walltime = clock.wallTime();
        time_unc = std::min(walltime, time_unc);
    }

    // compress state test
    double time_comp = std::numeric_limits<double>::max();
    M_performance.compress();
    for(std::size_t i = 0; i < n; ++i)
    {
        clock.start();
        res1 = std::move(M_performance * v);
        clock.stop();
        auto walltime = clock.wallTime();
        time_comp = std::min(walltime, time_comp);
    }

    // output
    std::cout << "----- Multiplication performance test results -----" << std::endl;
    std::cout << "Matrix read from file:     '" << filename_performance << "'" << std::endl;
    std::cout << "Storing method:            " << (
                ordering_performance == algebra::StorageOrder::ROW_WISE ? "'ROW-WISE'" : "'COLUMN-WISE'") << std::endl;
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Uncompressed best time:    " << time_unc << " microseconds" << std::endl;
    std::cout << "Compressed best time:      " << time_comp << " microseconds" << std::endl;
    std::cout << std::endl;


    //////// performing multiplication tests (matrix by vector but with Matrix-Matrix types) ///////////////////
    algebra::Matrix<var_type, ordering_mat> LM;
    algebra::Matrix<var_type, ordering_vector> RM;
    LM.parse_from_file(filename_mat);
    RM.parse_from_file(filename_vector);
    // compress matrices
    LM.compress();
    RM.compress();
    RM.uncompress();
    RM.compress();
    std::vector<var_type> M_res; // results
    M_res = std::move(LM * RM); // multiplication

    // output
    std::cout << "----- Multiplication Matrix-Matrix type test -----" << std::endl;
    std::cout << "Matrix read from file:     '" << filename_mat << "'" << std::endl;
    std::cout << "Vector read from file:     '" << filename_vector << "'" << std::endl;
    std::cout << "Matrix storing method:     " << (
                ordering_mat == algebra::StorageOrder::ROW_WISE ? "'ROW-WISE'" : "'COLUMN-WISE'") << std::endl;
    std::cout << "Vector storing method:     " << (
                ordering_vector == algebra::StorageOrder::ROW_WISE ? "'ROW-WISE'" : "'COLUMN-WISE'") << std::endl;
    // uncompress matrices to call the call operator that can print zeros (otherwise we get an error)                
    LM.uncompress();
    RM.uncompress();
    std::cout << std::endl;
    std::cout << "Matrix layout:                   " << std::endl;    
    for(std::size_t i = 0; i < LM.rows(); ++i)
    {
        std::cout << "                           | "; 
        for(std::size_t j = 0; j < LM.columns(); ++j)
        {
            std::cout << "(" << std::fixed << std::setprecision(1) << std::setw(4) << LM(i, j) << ") ";
        }
        std::cout << " |" << std::endl;
    }

    std::cout << "Vector layout:                   " << std::endl;    
    for(std::size_t i = 0; i < RM.rows(); ++i)
    {
        std::cout << "                           | ";
        std::cout << "(" << std::fixed << std::setprecision(1) << std::setw(4) << RM(i, 0) << ") |"; 
        std::cout << std::endl;
    }

    std::cout << "Resulting vector:" << std::endl;
    for(std::size_t i = 0; i < M_res.size(); ++i)
    {
        std::cout << "                           | ";
        std::cout << "(" << std::fixed << std::setprecision(1) << std::setw(4) << M_res[i] << ") |"; 
        std::cout << std::endl;
    }
    std::cout << std::endl;



    ///////// performing norm tests ////////////////////
    algebra::Matrix<var_type, ordering_norms> M_norms;
    M_norms.parse_from_file(filename_norms);

    std::cout << "--------------- Norms test results ----------------" << std::endl;
    std::cout << "Matrix read from file:     '" << filename_norms << "'" << std::endl;
    std::cout << "Storing method:            " << (
                ordering_norms == algebra::StorageOrder::ROW_WISE ? "'ROW-WISE'" : "'COLUMN-WISE'") << std::endl;
    // One norm
    M_norms.uncompress();
    var_type norm11 = M_norms.norm<norm_type1>();
    std::cout << "Norm type:                 " << (
                norm_type1 == algebra::NormType::One ? "'ONE'" : 
                    norm_type1 == algebra::NormType::Infinity ? "'INFINITY'" : "'FROBENIUS'") << std::endl;
    std::cout << "Uncompressed state" << std::endl;
    std::cout << "Norm:                      " << std::scientific << std::setprecision(16) << norm11 << std::endl;

    M_norms.compress();
    var_type norm12 = M_norms.norm<norm_type1>();
    std::cout << "Compressed state" << std::endl;
    std::cout << "Norm:                      " << std::scientific << std::setprecision(16) << norm12 << std::endl;
    std::cout << std::endl;
    // Infinity norm
    M_norms.uncompress();
    var_type norm21 = M_norms.norm<norm_type2>();
    std::cout << "Norm type:                 " << (
                norm_type2 == algebra::NormType::One ? "'ONE'" : 
                    norm_type2 == algebra::NormType::Infinity ? "'INFINITY'" : "'FROBENIUS'") << std::endl;
    std::cout << "Uncompressed state" << std::endl;
    std::cout << "Norm:                      " << std::scientific << std::setprecision(16) << norm21 << std::endl;

    M_norms.compress();
    var_type norm22 = M_norms.norm<norm_type2>();
    std::cout << "Compressed state" << std::endl;
    std::cout << "Norm:                      " << std::scientific << std::setprecision(16) << norm22 << std::endl;
    std::cout << std::endl;
    // Frobenius norm
    M_norms.uncompress();
    var_type norm31 = M_norms.norm<norm_type3>();
    std::cout << "Norm type:                 " << (
                norm_type3 == algebra::NormType::One ? "'ONE'" : 
                    norm_type3 == algebra::NormType::Infinity ? "'INFINITY'" : "'FROBENIUS'") << std::endl;
    std::cout << "Uncompressed state" << std::endl;
    std::cout << "Norm:                      " << std::scientific << std::setprecision(16) << norm31 << std::endl;

    M_norms.compress();
    var_type norm32 = M_norms.norm<norm_type3>();
    std::cout << "Compressed state" << std::endl;
    std::cout << "Norm:                      " << std::scientific << std::setprecision(16) << norm32 << std::endl;
    std::cout << std::endl;
    
    return 0;
}