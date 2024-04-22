#ifndef HH_DYNAMIC_MATRIX_HH
#define HH_DYNAMIC_MATRIX_HH

#include "MatrixTraits.hpp"
#include "MatrixBase.hpp"
#include <iostream>
#include<map>
#include<array>
namespace algebra
{
template<typename T>
class DynamicMatrix : public MatrixBase<T>
{

using elements_type = std::map<std::array<std::size_t, DIM>, T>; // just to ease notation

private:
    std::map<std::array<std::size_t, DIM>, T> elements{}; // elements

public:
    // leave to the compiler the default construtors
    DynamicMatrix() = default;

    DynamicMatrix(const std::size_t nrows, const std::size_t ncols)
    : MatrixBase<T>(nrows, ncols)
    {};

    // copy/move constructor
    DynamicMatrix(const std::size_t nrows, const std::size_t ncols, elements_type && elements_)
    : 
    MatrixBase<T>(nrows, ncols), 
    elements(std::forward<elements_type>(elements_)) 
    {};
    
    // call operator (non const version)
    T& operator()(const std::size_t i, const std::size_t j) override
    {
        if(i >= MatrixBase<T>::n_rows || j >= MatrixBase<T>::n_cols)
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
    const T operator()(const std::size_t i, const std::size_t j) const override
    {
        if(i >= MatrixBase<T>::n_rows || j >= MatrixBase<T>::n_cols)
            throw std::out_of_range("Invalid coordinates");

        auto iter = elements.find({i, j});
        
        if(iter == elements.end())
            return static_cast<T>(0.0);

        return iter->second;
    }

}; // end of class DynamicMatrix

} // end of namespace algebra


#endif