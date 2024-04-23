#ifndef HH_DYNAMIC_MATRIX_HH
#define HH_DYNAMIC_MATRIX_HH

#include "MatrixTraits.hpp"
#include <iostream>
#include<map>
#include<array>
namespace algebra
{
template<typename T>
class DynamicMatrix
{

using elements_type = std::map<std::array<std::size_t, DIM>, T>; // just to ease notation

public:
    elements_type elements{}; // elements

    ////// METHODS /////////
    // leave to the compiler the default construtors
    DynamicMatrix() = default;
    
    // copy/move constructor
    DynamicMatrix(elements_type && elements_)
    : 
    elements(std::forward<elements_type>(elements_)) 
    {};
    
    // call operator (non const version)
    T& operator()(const std::size_t i, const std::size_t j)
    {
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
        auto iter = elements.find({i, j});
        
        if(iter == elements.end())
            return static_cast<T>(0.0);

        return iter->second;
    }

    void resize(const std::size_t nrows, const std::size_t ncols)
    {
        elements_type new_elements;

        auto in_range = [nrows, ncols](const auto& element)
        {
            return element.first[0] < nrows && element.first[1] < ncols;
        };

        std::copy_if(elements.rbegin(), elements.rend(), std::inserter(new_elements, new_elements.end()), in_range);

        elements = std::move(new_elements);
    }

    // this method clears the stored matrix, it is needed to efficiently convert from dynamic to compress matrix
    void clear()
    {
        elements.clear();
    }

}; // end of class DynamicMatrix

} // end of namespace algebra


#endif