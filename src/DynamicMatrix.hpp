#ifndef HH_DYNAMIC_MATRIX_HH
#define HH_DYNAMIC_MATRIX_HH

#include "MatrixTraits.hpp"
#include <iostream>
#include<map>
#include<array>
namespace algebra
{
template<typename T, StorageOrder Order>
class DynamicMatrix
{

using elements_type = std::map<std::array<std::size_t, DIM>, T>; // just to ease notation

private:
    size_t n_rows = 0; // number of rows
    size_t n_cols = 0; // number of columns
    std::map<std::array<std::size_t, DIM>, T> elements; // elements

public:
    // leave to the compiler the default construtors
    DynamicMatrix() = default;

    // copy/move constructor
    DynamicMatrix(const std::size_t nrows, const std::size_t ncols, elements_type && elements_)
    : n_rows(nrows), n_cols(ncols), elements(std::forward<elements_type>(elements_)) 
    {};
    
    // call operator (non const version)
    T& operator()(const std::size_t i, const std::size_t j)
    {
        if(i > n_rows || j > n_cols)
            throw std::out_of_range("Invalid coordinates");

    
        auto iter = elements.find({i, j});

        if(iter == elements.end()) // allows to add a new element of type T
        {
            std::array<std::size_t, DIM> index{i, j};
            auto pair = elements.emplace(std::make_pair(index, T{})); // create a default element of type T at place (i, j)
            iter = pair.first; // pair returns a std::pair with first element the std::map and second 
                                // element a bool to accounts for possible problems during emplace() call
        }
        return iter->second;
    }

    // call operator (const version)
    const T operator()(const std::size_t i, const std::size_t j) const
    {
        if(i > n_rows || j > n_cols)
            throw std::out_of_range("Invalid coordinates");

        auto iter = elements.find({i, j});
        
        if(iter == elements.end())
            return static_cast<T>(0.0);

        return iter->second;
    }

    const std::size_t rows() const { return n_rows; }

    const std::size_t columns() const { return n_cols; }

}; // end of class DynamicMatrix

} // end of namespace algebra


#endif