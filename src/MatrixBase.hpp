#ifndef HH_MATRIX_BASE_HH
#define HH_MATRIX_BASE_HH

namespace algebra
{
    enum StorageOrder
    {
        ROW_WISE,
        COLUMN_WISE
    };

    template <typename T, StorageOrder Order>
    class Matrix
    {
    public:
        using value_type = T;

        Matrix(const std::size_t nrows, const std::size_t ncols)
            : n_rows(nrows), n_cols(ncols)
        {
        }

        Matrix(const std::size_t n)
            : Matrix(n, n)
        {
        }

        Matrix() = default;

        void resize(const std::size_t nrows, const std::size_t ncols)
        {
            n_rows = nrows;
            n_cols = ncols;

            if constexpr (Order == StorageOrder::ROW_WISE)
            {
                row_map.clear();
                row_map.reserve(nrows);
            }
            else
            {
                col_map.clear();
                col_map.reserve(ncols);
            }
        }

        const std::size_t rows() const { return n_rows; }
        const std::size_t cols() const { return n_cols; }

        const bool is_compressed() const
        {
            return (Order == StorageOrder::ROW_WISE ? row_map.empty() : col_map.empty());
        }

        void compress()
        {
            if (is_compressed())
                return;

            if constexpr (Order == StorageOrder::ROW_WISE)
            {
                std::vector<std::pair<std::array<std::size_t, 2>, T>> vec_map(row_map.begin(), row_map.end());
                row_map.clear();

                std::sort(vec_map.begin(), vec_map.end(),
                    [](const auto& lhs, const auto& rhs)
                    {
                        return lhs.first[0] < rhs.first[0] ||
                            (lhs.first[0] == rhs.first[0] && lhs.first[1] < rhs.first[1]);
                    });

                row_map.emplace_back(std::make_pair(vec_map[0].first, vec_map[0].second));

                for (std::size_t i = 1; i < vec_map.size(); ++i)
                {
                    if (row_map.back().first[0] != vec_map[i].first[0])
                    {
                        row_map.emplace_back(std::make_pair(vec_map[i].first, vec_map[i].second));
                    }
                    else
                    {
                        row_map.back().second += vec_map[i].second;
                    }
                }
            }
            else
            {
                std::vector<std::pair<std::array<std::size_t, 2>, T>> vec_map(col_map.begin(), col_map.end());
                col_map.clear();

                std::sort(vec_map.begin(), vec_map.end(),
                    [](const auto& lhs, const auto& rhs)
                    {
                        return lhs.first[1] < rhs.first[1] ||
                            (lhs.first[1] == rhs.first[1] && lhs.first[0] < rhs.first[0]);
                    });

                col_map.emplace_back(std::make_pair(vec_map[0].first, vec_map[0].second));

                for (std::size_t i = 1; i < vec_map.size(); ++i)
                {
                    if (col_map.back().first[1] != vec_map[i].first[1])
                    {
                        col_map.emplace_back(std::make_pair(vec_map[i].first, vec_map[i].second));
                    }
                    else
                    {
                        col_map.back().second += vec_map[i].second;
                    }
                }
            }
        }

        void uncompress()
        {
            if (!is_compressed())
                return;

            if constexpr (Order == StorageOrder::ROW_WISE)
            {
                row_map.reserve(n_rows * n_cols);

                for (std::size_t i = 0; i < n_rows; ++i)
                {
                    for (std::size_t j = 0; j < n_cols; ++j)
                    {
                        row_map.emplace_back(std::make_pair(std::array<std::size_t, 2>{i, j}, T{}));
                    }
                }
            }
            else
            {
                col_map.reserve(n_rows * n_cols);

                for (std::size_t i = 0; i < n_rows; ++i)
                {
                    for (std::size_t j = 0; j < n_cols; ++j)
                    {
                        col_map.emplace_back(std::make_pair(std::array<std::size_t, 2>{i, j}, T{}));
                    }
                }
            }
        }

        T& operator()(const std::size_t i, const std::size_t j)
        {
            if (!(i < n_rows && j < n_cols))
                throw std::out_of_range("Matrix index out of range");

            if constexpr (Order == StorageOrder::ROW_WISE)
            {
                if (is_compressed())
                {
                    auto iter = row_map.find(std::array<std::size_t, 2>{i, j});

                    if (iter == row_map.end())
                    {
                        return row_map.emplace(std::make_pair(std::array<std::size_t, 2>{i, j}, T{})).first->second;
                    }

                    return iter->second;
                }
                else
                {
                    return row_map[i * n_cols + j];
                }
            }
            else
            {
                if (is_compressed())
                {
                    auto iter = col_map.find(std::array<std::size_t, 2>{i, j});

                    if (iter == col_map.end())
                    {
                        return col_map.emplace(std::make_pair(std::array<std::size_t, 2>{i, j}, T{})).first->second;
                    }

                    return iter->second;
                }
                else
                {
                    return col_map[j * n_rows + i];
                }
            }
        }

        const T operator()(const std::size_t i, const std::size_t j) const
        {
            if (!(i < n_rows && j < n_cols))
                return T{};

            if constexpr (Order == StorageOrder::ROW_WISE)
            {
                if (is_compressed())
                {
                    auto iter = row_map.find(std::array<std::size_t, 2>{i, j});

                    if (iter == row_map.end())
                        return T{};

                    return iter->second;
                }
                else
                {
                    return row_map[i * n_cols + j];
                }
            }
            else
            {
                if (is_compressed())
                {
                    auto iter = col_map.find(std::array<std::size_t, 2>{i, j});

                    if (iter == col_map.end())
                        return T{};

                    return iter->second;
                }
                else
                {
                    return col_map[j * n_rows + i];
                }
            }
        }

        template <typename U>
        friend Matrix<std::common_type_t<T, U>, Order> operator*(const Matrix<T, Order>& A, const std::vector<U>& v)
        {
            if (A.cols() != v.size())
                throw std::invalid_argument("Matrix-vector multiplication: invalid dimensions");

            Matrix<std::common_type_t<T, U>, Order> C(A.rows(), 1);
        }

    private:
        std::size_t n_rows;
        std::size_t n_cols;
        bool is_compressed;
        CompressedMatrix<T, Order> ComprMat;
        UncompressedMatrix<T, Order> UncomprMat;
    }; // end of class Matrix
} // end namespace algebra
#endif