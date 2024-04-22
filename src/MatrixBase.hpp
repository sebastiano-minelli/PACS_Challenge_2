#ifndef HH_MATRIX_BASE_HH
#define HH_MATRIX_BASE_HH

#include<iostream>
namespace algebra
{

template<typename T>
class MatrixBase
{
protected:
    std::size_t n_rows = 0;
    std::size_t n_cols = 0;

public:

    MatrixBase() = default;

    MatrixBase(const std::size_t nrows, const std::size_t ncols)
        : n_rows(nrows), n_cols(ncols)
    {};

    virtual T& operator()(const std::size_t i, const std::size_t j) = 0; // virtual call operator

    virtual const T operator()(const std::size_t i, const std::size_t j) const = 0; // virtual const call operator

    // inline const std::size_t rows() const { return n_rows; }

    // inline const std::size_t columns() const { return n_cols; }

}; // end of class MatrixBase

} // end namespace algebra

#endif