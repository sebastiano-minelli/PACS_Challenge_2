#ifndef HH_MATRIX_TRAITS_HH
#define HH_MATRIX_TRAITS_HH

#include<iostream>
namespace algebra
{
// Storage order enumerator
enum StorageOrder
{
    ROW_WISE,
    COLUMN_WISE
};

inline constexpr std::size_t DIM = 2; // Dimension of the tuple that describes a matrix element

}

#endif