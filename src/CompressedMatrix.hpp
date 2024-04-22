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
    T operator()(const std::size_t i, const std::size_t j)
    {
        if constexpr (Order == StorageOrder::ROW_WISE)
        {
            if (i < 0 || i >= n_rows || j < 0 || j >= n_cols) 
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
            if (i < 0 || i >= n_rows || j < 0 || j >= n_cols) 
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

    // call operator (const version)
    const T operator()(const std::size_t i, const std::size_t j) const
    {
        if constexpr (Order == StorageOrder::ROW_WISE)
        {
            if (i < 0 || i >= n_rows || j < 0 || j >= n_cols) 
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
            if (i < 0 || i >= n_rows || j < 0 || j >= n_cols) 
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

    /*
    void set(std::size_t i, std::size_t j, T value)
    {
        if constexpr (Order == StorageOrder::ROW_WISE)
        {
            auto outer = get_outer_index(i, j);
            auto inner_begin = inner_indexes.begin() + i;
            auto inner_end = inner_indexes.begin() + i + 1;

            auto iter = std::lower_bound(inner_begin, inner_end, outer);

            if (iter != inner_end && *iter == outer)
            {
                values[std::distance(inner_indexes.begin(), iter)] += value;
            }
            else
            {
                inner_indexes.insert(iter, outer);
                outer_indexes.push_back(j);
                values.push_back(value);
            }
        }
        else
        {
            auto inner = get_inner_index(i, j);
            auto outer_begin = outer_indexes.begin() + inner;
            auto outer_end = outer_indexes.begin() + inner + 1;

            auto iter = std::lower_bound(outer_begin, outer_end, i);

            if (iter != outer_end && *iter == i)
            {
                values[std::distance(outer_indexes.begin(), iter)] += value;
            }
            else
            {
                outer_indexes.insert(iter, i);
                inner_indexes[j]++;
                values.insert(values.begin() + inner, value);
            }
        }
    }
    */

};

} // end of namespace algebra



#endif // HH_COMPRESSED_MATRIX_HH