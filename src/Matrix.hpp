#ifndef HH_MATRIX_BASE_HH
#define HH_MATRIX_BASE_HH

#include<memory>
#include "DynamicMatrix.hpp"
#include "CompressedMatrix.hpp"
namespace algebra
{
template <typename T, StorageOrder Order>
class Matrix
{
private:
    bool is_compressed = false; // true = comrpessed storing
    std::size_t n_rows = static_cast<std::size_t>(0); // number of rows
    std::size_t n_cols = static_cast<std::size_t>(0); // number of columns
    DynamicMatrix dynamic_mat{}; // dynamic matrix data
    CompressedMatrix compressed_mat{}; // compressed matrix data


public:
    // create a matrix, set it to uncompressed state (even if it is empty)
    Matrix(const std::size_t nrows, const std::size_t ncols)
        : n_rows(nrows), n_cols(ncols)
    {
        DynamicMatrix{};
        CompressedMatrix{};
        is_compressed = false;
    };

    // create a matrix, set it to uncompressed state
    using elements_type = std::map<std::array<std::size_t, DIM>, T>; // just to ease notation
    Matrix(const std::size_t nrows, const std::size_t ncols, elements_type && elements_)
    : 
    is_compressed(false),
    n_rows(nrows), 
    n_cols(ncols), 
    dynamic_mat(std::move(elements_))
    {
        CompressedMatrix{};
    };

    void resize(const std::size_t nrows, const std::size_t ncols)
    {
        if(!is_compressed && (nrows < n_rows || ncols < n_cols)) // if where have a dynamic matrix we need to call the resize() method
        {
            n_rows = nrows;
            n_cols = ncols;
            dynamic_mat.resize(n_rows, n_cols);
        }
        else if(!is_compressed) // if it is dynamic but we enlarge it no need to take out elements
        {
            n_rows = nrows;
            n_cols = ncols;
        }
        else // if(is_compressed)
        {
            std::cout << "Attention: the matrix will be leaved in an uncompressed state" << std::endl;
            n_rows = nrows;
            n_cols = ncols;
            this->uncompress(); // uncompress takes care of adding or deleting elements
        }
    }

    const std::size_t rows() const { return n_rows; }

    const std::size_t columns() const { return n_cols; }


    const bool is_compressed() const
    {
        return is_compressed;
    }

    void compress()
    {
        if(is_compressed)
            return;
        //else

        // Clear the current compressed matrix data (just to be shure)
        compressed_mat.clear();

        compressed_mat.values.reserve(dynamic_mat.size()); // reserving space

        const auto first_row = dynamic_mat.elements.cbegin()->first[0];
        const auto last_row = dynamic_mat.elements.crbegin()->first[0];
        const auto first_col = dynamic_mat.elements.cbegin()->first[1];
        const auto last_col = dynamic_mat.elements.crbegin()->first[1];

        if constexpr (Order == StorageOrder::ROW_WISE)
        {
            // reserve space
            compressed_mat.inner_indexes.reserve(last_row - first_row + 2); 
            compressed_mat.outer_indexes.reserve(dynamic_mat.size());

            // store all column indexes inside outer_indexes
            std::ranges::transform(dynamic_mat.elements.cbegin(),
                                dynamic_mat.elements.cend(), 
                                std::back_inserter(compressed_mat.outer_indexes), 
                                [](auto& pair) { return pair->first[1]; });

            // store row indexes inside inner_indexes
            compressed_mat.inner_indexes.push_back(0); // first element is always 0
            
            for(std::size_t i = first_row; i <= last_row; ++i)
            {
                auto low_bound = dynamic_mat.elements.lower_bound({i, 0});
                auto up_bound = dynamic_mat.elements.upper_bound({i + 1, 0});
                auto distance = std::ranges::distance(low_bound, up_bound); // computing how many elements there are in a row
                compressed_mat.inner_indexes.push_back(distance);
            }           

        }
        else // if Order == StorageOrder::COLUMN_WISE
        {
            // reserve space
            compressed_mat.inner_indexes.reserve(last_col - first_col + 2); 
            compressed_mat.outer_indexes.reserve(dynamic_mat.size());

            // store all row indexes inside outer_indexes
            std::ranges::transform(dynamic_mat.elements.cbegin(),
                                dynamic_mat.elements.cend(), 
                                std::back_inserter(compressed_mat.outer_indexes), 
                                [](auto& pair) { return pair->first[0]; });

            // store column indexes inside inner_indexes
            compressed_mat.inner_indexes.push_back(0); // first element is always 0
            
            for(std::size_t i = first_col; i <= last_col; ++i)
            {
                auto low_bound = dynamic_mat.elements.lower_bound({0, i});
                auto up_bound = dynamic_mat.elements.upper_bound({0, i + 1});
                auto distance = std::ranges::distance(low_bound, up_bound); // computing how many elements there are in a row
                compressed_mat.inner_indexes.push_back(distance);
            }

        }

        // move elements from the map to the compressed matrix
            std::ranges::transform(dynamic_mat.elements.begin(),
                                dynamic_mat.elements.end(), 
                                std::back_inserter(compressed_mat.values), 
                                [](auto& pair) { return std::move(pair->second); });        

        // Change compressed state
        is_compressed = true;

        // Clear the dynamic matrix data
        dynamic_mat.clear();
    }
    

    void uncompress()
    {
        if (!is_compressed())
            return;
        // else

        // Clear the current dynamic matrix data (just to be shure)
        dynamic_mat.clear();

        const auto outer_it = compressed_mat.outer_indexes.cbegin();
        const auto value_it = compressed_mat.values.cbegin();

        if constexpr (Order == StorageOrder::ROW_WISE)
        {
            for(std::size_t i = 0; i < compressed_mat.inner_indexes.size(); ++i)
            {
                for(std::size_t j = 0; j < compressed_mat.inner_indexes[i]; ++j)
                {
                    dynamic_mat.insert(std::make_pair(std::array<std::size_t, DIM>{i, *outer_it}, *value_it));
                    ++outer_it;
                    ++value_it;
                }
            }
        }
        else // if Order == StorageOrder::COLUMN_WISE
        {
            for(std::size_t j = 0; j < compressed_mat.inner_indexes.size(); ++j)
            {
                for(std::size_t i = 0; i < compressed_mat.inner_indexes[j]; ++i)
                {
                    dynamic_mat.insert(std::make_pair(std::array<std::size_t, DIM>{*outer_it, j}, *value_it));
                    ++outer_it;
                    ++value_it;
                }
            }
        }
    }

    T& operator()(const std::size_t i, const std::size_t j)
    {
        if (i > n_rows || j > n_cols)
            throw std::out_of_range("Matrix index out of range");
        if(is_compressed)
            return compressed_mat(i, j);
        return dynamic_mat(i, j);
    }

    const T operator()(const std::size_t i, const std::size_t j) const
    {
        if (i > n_rows || j > n_cols)
            throw std::out_of_range("Matrix index out of range");
        if(is_compressed)
            return const compressed_mat(i, j);
        return const dynamic_mat(i, j);
    }

    template <typename U>
    friend Matrix<std::common_type_t<T, U>, Order> operator*(const Matrix<T, Order>& A, const std::vector<U>& v)
    {
        if (A.cols() != v.size())
            throw std::invalid_argument("Matrix-vector multiplication: invalid dimensions");

        Matrix<std::common_type_t<T, U>, Order> C(A.rows(), 1);
    }

}; // end of class Matrix

} // end namespace algebra

#endif