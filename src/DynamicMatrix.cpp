#ifndef HH_DYNAMIC_MATRIX_HH
#define HH_DYNAMIC_MATRIX_HH

namespace algebra
{
template<typename T, StorageOrder Order = StorageOrder::ROW_WISE>
class Matrix
{
public:
    Matrix(size_t nrows = 0, size_t ncols = 0) : n_rows(nrows), n_cols(ncols) {}

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
    friend std::ostream& operator<<(std::ostream& os, const Matrix<U, Order>& A);

private:
    size_t n_rows;
    size_t n_cols;
    std::map<std::array<std::size_t, 2>, T, std::less<std::array<std::size_t, 2>>> map;
};

template<typename T, StorageOrder Order>
std::ostream& operator<<(std::ostream& os, const Matrix<T, Order>& A)
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