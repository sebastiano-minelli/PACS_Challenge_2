#ifndef HH_MATRIX_HH
#define HH_MATRIX_HH

namespace algebra
{

enum class StorageOrder
{
    ROW_WISE = 0,
    COL_WISE = 1
};

template<class T, StorageOrder Order = StorageOrder::ROW_WISE>
class Matrix
{
public:
    Matrix<T, StorageOrder>(size_t nrows, size_t ncols) :
    n_rows(nrows), n_cols(ncols)
    {};

private:
    size_t n_rows; // number of rows
    size_t n_cols; // number of columns


};

} // end of namespace algebra

#endif // HH_MATRIX_HH