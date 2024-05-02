#ifndef HH_MATRIX_UTILITIES_HH
#define HH_MATRIX_UTILITIES_HH

#include <iostream>
#include <map>
namespace algebra
{
// Storage order enumerator
enum StorageOrder
{
    ROW_WISE,
    COLUMN_WISE
};

// Enumerator that lists the possible norms
enum NormType
{
    One,
    Infinity,
    Frobenius
};

inline constexpr std::size_t DIM = 2; // Dimension of the tuple that describes a matrix element

/// Type specialization for DynamicMatrix class elements ///
// we need to change the standard comparison operator for arrays if we traverse the std::map COLUMN_WISE to have the nice features 
// of std::map, provide a struct that selects the right elements_type automatically
template<typename T, StorageOrder O>
struct ElementsTypeSelection;

// specialization for StorageOrder::ROW_WISE
template<typename T>
struct ElementsTypeSelection<T, StorageOrder::ROW_WISE>
{
    using type = std::map<std::array<std::size_t, DIM>, T>;
};

// specialization for StorageOrder::COLUMN_WISE
template<typename T>
struct ElementsTypeSelection<T, StorageOrder::COLUMN_WISE>
{
    // define comparison operator
    struct compare
    {
        bool operator()(const std::array<std::size_t, DIM>& a, const std::array<std::size_t, DIM>& b) const
        {
            for (std::size_t i = DIM; i > 0; --i) 
            {
                if (a[i - 1] < b[i - 1])
                    return true;
                if (a[i - 1] > b[i - 1])
                    return false;
            }
            return false; // arrays are equal
        }; // notice that here I used DIM that is actually known, I could have used directly DIM = 2 but I think (hope) that with
        // optimizations activated the code will be equivalent
    };

    using type = std::map<std::array<std::size_t, DIM>, T, compare>;
};


}

#endif