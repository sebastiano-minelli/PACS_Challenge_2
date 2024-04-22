#ifndef HH_DYNAMIC_MATRIX_HH
#define HH_DYNAMIC_MATRIX_HH

#include <iostream>
namespace algebra
{
// Storage order enumerator
enum StorageOrder
{
    ROW_WISE,
    COLUMN_WISE
};

inline constexpr std::size_t DIM = 2; // Dimension of the tuple that describes a matrix element

template<typename T, StorageOrder Order = StorageOrder::ROW_WISE>
class DynamicMatrix
{

private:
    size_t n_rows; // numer of rows
    size_t n_cols; // number of columns
    std::map<std::array<std::size_T, DIM>, T> elements; // elements

public:
    // leave to the compiler the default construtors
    
    T& operator()(const std::size_t i, const std::size_t j)
    {
        if (!(i < n_rows && j < n_cols))
            throw std::out_of_range("Matrix access: invalid coordinates");

        if constexpr (Order == StorageOrder::ROW_WISE)
        {
            auto iter = map.find(std::array<std::size_t, 2>{i, j});

            if (iter == map.end())
                return map.emplace(std::make_pair(std::array<std::size_t, 2>{i, j}, T{})).first->second;

            return iter->second;
        }
        else
        {
            auto iter = map.find(std::array<std::size_t, 2>{j, i});

            if (iter == map.end())
                return map.emplace(std::make_pair(std::array<std::size_t, 2>{j, i}, T{})).first->second;

            return iter->second;
        }
    }

    const T& operator()(const std::size_t i, const std::size_t j) const
    {
        if (!(i < n_rows && j < n_cols))
            throw std::out_of_range("Matrix access: invalid coordinates");

        if constexpr (Order == StorageOrder::ROW_WISE)
        {
            auto iter = map.find(std::array<std::size_t, 2>{i, j});

            if (iter == map.end())
                return T{};

            return iter->second;
        }
        else
        {
            auto iter = map.find(std::array<std::size_t, 2>{j, i});

            if (iter == map.end())
                return T{};

            return iter->second;
        }
    }

    size_t rows() const { return n_rows; }

    size_t cols() const { return n_cols; }

    template<typename U>
    friend std::ostream& operator<<(std::ostream& os, const DynamicMatrix<U, Order>& A);

}; // end of class DynamicMatrix

template<typename T, StorageOrder Order>
std::ostream& operator<<(std::ostream& os, const DynamicMatrix<T, Order>& A)
{
    for (size_t i = 0; i < A.rows(); ++i)
    {
        for (size_t j = 0; j < A.cols(); ++j)
            os << A(i, j) << "\t";

        os << std::endl;
    }

    return os;
}

} // end of namespace algebra


#endif