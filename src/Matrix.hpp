#ifndef HH_MATRIX_BASE_HH
#define HH_MATRIX_BASE_HH

#include <fstream>
#include <sstream>
#include <string>
#include <memory>
#include "DynamicMatrix.hpp"
#include "CompressedMatrix.hpp"
namespace algebra
{
template <typename T, StorageOrder Order>
class Matrix
{
private:
    bool compressed = false; // true = compressed storing
    std::size_t n_rows = static_cast<std::size_t>(0); // number of rows
    std::size_t n_cols = static_cast<std::size_t>(0); // number of columns
    DynamicMatrix<T, Order> dynamic_mat{}; // dynamic matrix data
    CompressedMatrix<T, Order> compressed_mat{}; // compressed matrix data


public:

    Matrix() = default;

    // create a matrix, set it to uncompressed state (even if it is empty)
    Matrix(const std::size_t nrows, const std::size_t ncols)
    : 
    n_rows(nrows), n_cols(ncols)
    {};

    // create a matrix, set it to uncompressed state
    using elements_type = typename algebra::ElementsTypeSelection<T, Order>::type; // selecting the right type
    Matrix(const std::size_t nrows, const std::size_t ncols, elements_type && elements_)
    : 
    compressed(false),
    n_rows(nrows), 
    n_cols(ncols), 
    dynamic_mat(std::move(elements_))
    {};

    void resize(const std::size_t nrows, const std::size_t ncols)
    {
        if(!compressed && (nrows < n_rows || ncols < n_cols)) // if where have a dynamic matrix we need to call the resize() method
        {
            n_rows = nrows;
            n_cols = ncols;
            dynamic_mat.resize(n_rows, n_cols);
        }
        else if(!compressed) // if it is dynamic but we enlarge it no need to take out elements
        {
            n_rows = nrows;
            n_cols = ncols;
        }
        else // if(compressed)
        {
            std::cout << "Attention: the matrix will be leaved in an uncompressed state" << std::endl;
            n_rows = nrows;
            n_cols = ncols;
            this->uncompress();
            dynamic_mat.resize(n_rows, n_cols);
        }
    }

    const std::size_t rows() const { return n_rows; }

    const std::size_t columns() const { return n_cols; }


    const bool is_compressed() const
    {
        return compressed;
    }

    void compress()
    {
        if(compressed)
            return;
        //else

        // Clear the current compressed matrix data (just to be shure)
        compressed_mat.clear();

        compressed_mat.values.reserve(dynamic_mat.elements.size()); // reserving space

        if constexpr (Order == StorageOrder::ROW_WISE)
        {
            // reserve space
            compressed_mat.inner_indexes.reserve(n_rows + 1);
            compressed_mat.outer_indexes.reserve(dynamic_mat.elements.size());

            // store all column indexes inside outer_indexes
            std::transform(dynamic_mat.elements.cbegin(),
                                dynamic_mat.elements.cend(), 
                                std::back_inserter(compressed_mat.outer_indexes), 
                                [](auto& pair) { return pair.first[1]; });

            // store row indexes inside inner_indexes
            compressed_mat.inner_indexes.push_back(0); // first element is always 0
            
            std::size_t distance = static_cast<std::size_t>(0);
            std::cout << "distance = " << distance << std::endl;
            for(std::size_t i = 0; i < n_rows; ++i)
            {
                auto low_bound = dynamic_mat.elements.lower_bound({i, 0});
                auto up_bound = dynamic_mat.elements.upper_bound({i, n_cols});
                std::cout << "distance ranges: " << std::ranges::distance(low_bound, up_bound) << std::endl;
                distance += std::ranges::distance(low_bound, up_bound); // computing how many elements there are in a row
                std::cout << "distance = " << distance << std::endl;
                compressed_mat.inner_indexes.push_back(distance);
            }     

        }
        else // if Order == StorageOrder::COLUMN_WISE
        {
            // reserve space
            compressed_mat.inner_indexes.reserve(n_cols + 1); 
            compressed_mat.outer_indexes.reserve(dynamic_mat.elements.size());

            // store all row indexes inside outer_indexes
            std::transform(dynamic_mat.elements.cbegin(),
                                dynamic_mat.elements.cend(), 
                                std::back_inserter(compressed_mat.outer_indexes), 
                                [](auto& pair) { return pair.first[0]; });

            // store column indexes inside inner_indexes
            compressed_mat.inner_indexes.push_back(0); // first element is always 0
            std::size_t distance = static_cast<std::size_t>(0);
            for(std::size_t i = 0; i < n_cols; ++i)
            {
                auto low_bound = dynamic_mat.elements.lower_bound({0, i});
                auto up_bound = dynamic_mat.elements.upper_bound({n_rows, i});
                distance += std::ranges::distance(low_bound, up_bound); // computing how many elements there are in a row
                compressed_mat.inner_indexes.push_back(distance);
            }

        }

        // move elements from the map to the compressed matrix
        std::transform(dynamic_mat.elements.begin(),
                            dynamic_mat.elements.end(), 
                            std::back_inserter(compressed_mat.values), 
                            [](auto& pair) { return std::move(pair.second); });      

        // Change compressed state
        compressed = true;

        // Clear the dynamic matrix data
        dynamic_mat.clear();
    }
    

    void uncompress()
    {
        if (!compressed)
            return;
        // else

        // Clear the current dynamic matrix data (just to be shure)
        dynamic_mat.clear();

        auto outer_it = compressed_mat.outer_indexes.cbegin();
        auto value_it = compressed_mat.values.cbegin();

        if constexpr (Order == StorageOrder::ROW_WISE)
        {
            for(std::size_t i = 0; i < compressed_mat.inner_indexes.size(); ++i)
            {
                for(std::size_t j = 0; j < compressed_mat.inner_indexes[i + 1] - compressed_mat.inner_indexes[i]; ++j)
                {
                    dynamic_mat.elements.insert(std::make_pair(std::array<std::size_t, DIM>{i, *outer_it}, *value_it));
                    ++outer_it;
                    ++value_it;
                }
            }
        }
        else // if Order == StorageOrder::COLUMN_WISE
        {
            for(std::size_t j = 0; j < compressed_mat.inner_indexes.size(); ++j)
            {
                for(std::size_t i = 0; i < compressed_mat.inner_indexes[j + 1] - compressed_mat.inner_indexes[j]; ++i)
                {
                    dynamic_mat.elements.insert(std::make_pair(std::array<std::size_t, DIM>{*outer_it, j}, *value_it));
                    ++outer_it;
                    ++value_it;
                }
            }
        }

        // Change compressed state
        compressed = false;

        // Clear the compress matrix data
        compressed_mat.clear();
    }

    T& operator()(const std::size_t i, const std::size_t j)
    {
        if (i > n_rows || j > n_cols)
            throw std::out_of_range("Matrix index out of range");
        if(compressed)
            return compressed_mat(i, j);
        return dynamic_mat(i, j);
    }

    const T operator()(const std::size_t i, const std::size_t j) const
    {
        if (i > n_rows || j > n_cols)
            throw std::out_of_range("Matrix index out of range");
        if(compressed)
            return compressed_mat(i, j);
        return dynamic_mat(i, j);
    }

    friend std::vector<T> operator*(Matrix<T, Order>& mat, std::vector<T> & v)
    {
        if(mat.n_cols != v.size())
            throw std::invalid_argument("Matrix-vector multiplication: invalid dimensions");
        //else
        const auto dim = mat.n_rows;
        std::vector<T> result(dim); // initialized to T{}

        if(mat.compressed)
        {
            auto value_it = mat.compressed_mat.values.cbegin();
            const auto inner_index_dim = mat.compressed_mat.inner_indexes.size();
            auto outer_it = mat.compressed_mat.outer_indexes.cbegin();

            if constexpr (Order == StorageOrder::ROW_WISE)
            {                           
                for(std::size_t i = 0; i < inner_index_dim - 1; ++i)
                {
                    for(std::size_t j = 0; j < mat.compressed_mat.inner_indexes[i + 1] - mat.compressed_mat.inner_indexes[i]; ++j)
                    {
                        result[i] += (*value_it) * v[(*outer_it)];
                        ++value_it;
                        ++outer_it;
                    }
                }
            }
            else // if mat.Order == StorageOrder::COLUMN_WISE
            {
                for(std::size_t i = 0; i < inner_index_dim; ++i)
                {
                    for(std::size_t j = 0; j < mat.compressed_mat.inner_indexes[i + 1] - mat.compressed_mat.inner_indexes[i]; ++j)
                    {    
                        result[*outer_it] += (*value_it) * v[i];
                        ++value_it;
                        ++outer_it;
                    }
                }
            }
        }
        else // if mat.compressed = false
        {
            for(auto & element : mat.dynamic_mat.elements)
                result[element.first[0]] += element.second * v[element.first[1]];                
        }

        return result;
    }

    void parse_from_file(const std::string & filename)
    {
        // Check that the file exists
        std::ifstream file_check(filename);
        if(!file_check.is_open())
            throw std::runtime_error("Error: file not found");

        n_rows = 0;
        n_cols = 0;
        dynamic_mat.clear();
        compressed_mat.clear();
        compressed = false;

        // Parse as a dynamic matrix
        
        elements_type file_mat;

        std::ifstream file(filename);
        std::string line;

        // Skip the first line (that is just a comment)
        std::getline(file, line);

        // Read matrix dimensions
        if(std::getline(file, line)) 
        {
            std::istringstream word(line);
            word >> n_rows; 
            word >> n_cols;
            // I skip the number of non-zero elements since I don't need it to store the dynamic matrix
        }

        // Read the matrix data
        while(std::getline(file, line) )
        {
            // If the line starts with '%' we go on
            if(!line.empty() && line[0] == '%') // unnecessary probably since there are no comments
                continue;
            std::istringstream word(line);
            std::size_t i;
            std::size_t j;
            T value;
            if (!(word >> i >> j >> value))
                throw std::runtime_error("Error while reading the file");

            // Insert the element in the map
            std::cout << "element (" << i << ", " << j << ") = " << value << std::endl;
            file_mat.insert(std::make_pair(std::array<std::size_t, DIM>{i - 1, j - 1}, value));
        }

        dynamic_mat = DynamicMatrix<T, Order>(std::move(file_mat));
        std::cout << "Number of elements: " << dynamic_mat.elements.size() << std::endl;
    }

}; // end of class Matrix

} // end namespace algebra

#endif