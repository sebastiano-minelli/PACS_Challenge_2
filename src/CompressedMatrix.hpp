#ifndef HH_COMPRESSED_MATRIX_HH
#define HH_COMPRESSED_MATRIX_HH

#include "MatrixTraits.hpp"
#include <iostream>
#include <vector>

namespace algebra 
{
template<typename T, StorageOrder Order>
class CompressedMatrix
{
private:
    std::size_t n_rows = 0; // number of rows
    std::size_t n_cols = 0; // number of columns

    std::vector<std::size_t> inner_indexes; // vector that stores the inned indexes
    std::vector<std::size_t> outer_indexes; // vector that stores the outer indexes
    std::vector<T> values; // vector that stores the non zero values

public:

    CompressedMatrix() = default;

    CompressedMatrix(const std::size_t nrows, 
                    const std::size_t ncols,
                    std::vector<std::size_t> && inner_indexes_,
                    std::vector<std::size_t> && outer_indexes_, 
                    std::vector<T> && values_) 
    :
    n_rows(nrows), 
    n_cols(ncols),
    inner_indexes(std::forward<std::vector<std::size_t>>(inner_indexes_)),
    outer_indexes(std::forward<std::vector<std::size_t>>(outer_indexes_)),
    values(std::forward<std::vector<T>>(values_))
    {};

    // call operator (non const version)
    T& operator()(const std::size_t i, const std::size_t j)
    {
        if constexpr (Order == StorageOrder::ROW_WISE)
        {
            if (i >= n_rows || j >= n_cols) 
                throw std::out_of_range("Invalid coordinates");
            else
            {
                // compute range of elements to look for
                const std::size_t row_start = inner_indexes[i]; 
                const std::size_t row_end = inner_indexes[i + 1];

                // Look for the column index j within the range
                for (std::size_t k = row_start; k < row_end; ++k) 
                    if(outer_indexes[k] == j)
                        return values[k];
                throw std::out_of_range("Cannot assign to an invalid index");
            }
        }
        else // if Order == StorageOrder::COLUMN_WISE
        {
            if (i >= n_rows || j >= n_cols) 
                throw std::out_of_range("Invalid coordinates");
            else
            {
                // compute range of elements to look for
                const std::size_t col_start = inner_indexes[j]; 
                const std::size_t col_end = inner_indexes[j + 1];

                // Look for the column index j within the range
                for (std::size_t k = col_start; k < col_end; ++k) 
                    if(outer_indexes[k] == i)
                        return values[k];
                throw std::out_of_range("Cannot assign to an invalid index");
            }
        }
    }

    // call operator (const version)
    const T operator()(const std::size_t i, const std::size_t j) const
    {
        if constexpr (Order == StorageOrder::ROW_WISE)
        {
            if (i >= n_rows || j >= n_cols) 
                throw std::out_of_range("Invalid coordinates");
            else
            {
                // compute range of elements to look for
                const std::size_t row_start = inner_indexes[i]; 
                const std::size_t row_end = inner_indexes[i + 1];

                // Look for the column index j within the range
                for (std::size_t k = row_start; k < row_end; ++k) 
                    if(outer_indexes[k] == j)
                        return values[k];
                return static_cast<T>(0.0);
            }
        }
        else // if Order == StorageOrder::COLUMN_WISE
        {
            if (i >= n_rows || j >= n_cols) 
                throw std::out_of_range("Invalid coordinates");
            else
            {
                // compute range of elements to look for
                const std::size_t col_start = inner_indexes[j]; 
                const std::size_t col_end = inner_indexes[j + 1];

                // Look for the column index j within the range
                for (std::size_t k = col_start; k < col_end; ++k) 
                    if(outer_indexes[k] == i)
                        return values[k];
                return static_cast<T>(0.0);
            }
        }
    }

    const std::size_t rows() const { return n_rows; }

    const std::size_t columns() const { return n_cols; }    

};

} // end of namespace algebra



#endif // HH_COMPRESSED_MATRIX_HH