#ifndef HH_COMPRESSED_MATRIX_HH
#define HH_COMPRESSED_MATRIX_HH


template<typename T, StorageOrder Order>
class CompressedMatrix
{
private:
    //
    std::size_t n_rows;
    std::size_t n_cols;

    std::vector<std::size_t> inner_indexes;
    std::vector<std::size_t> outer_indexes;
    std::vector<T> values;


    std::size_t get_outer_index(std::size_t i, std::size_t j) const
    {
        if (i == 0)
            return 0;

        return inner_indexes[i - 1] + j - (inner_indexes[i - 1] == 0 ? 0 : outer_indexes[inner_indexes[i - 1] - 1] + 1);
    }

    std::size_t get_inner_index(std::size_t i, std::size_t j) const
    {
        return outer_indexes[inner_indexes[j] + (i > outer_indexes[inner_indexes[j]] ? 1 : 0) - (j > 0 && outer_indexes[inner_indexes[j - 1]] == outer_indexes[inner_indexes[j]] ? inner_indexes[j - 1] : 0)];
    }

public:

    public:
    CompressedMatrix(std::size_t nrows, std::size_t ncols) :
        n_rows(nrows), n_cols(ncols),
        inner_indexes(nrows + 1, 0), outer_indexes(), values()
    {}

    void set(std::size_t i, std::size_t j, T value)
    {
        if (value == T{})
            return;

        if constexpr (Order == StorageOrder::ROW_WISE)
        {
            auto outer = get_outer_index(i, j);
            auto inner_begin = inner_indexes.begin() + i;
            auto inner_end = inner_indexes.begin() + i + 1;

            auto iter = std::lower_bound(inner_begin, inner_end, outer);

            if (iter != inner_end && *iter == outer)
            {
                values[std::distance(inner_indexes.begin(), iter)] += value;
            }
            else
            {
                inner_indexes.insert(iter, outer);
                outer_indexes.push_back(j);
                values.push_back(value);
            }
        }
        else
        {
            auto inner = get_inner_index(i, j);
            auto outer_begin = outer_indexes.begin() + inner;
            auto outer_end = outer_indexes.begin() + inner + 1;

            auto iter = std::lower_bound(outer_begin, outer_end, i);

            if (iter != outer_end && *iter == i)
            {
                values[std::distance(outer_indexes.begin(), iter)] += value;
            }
            else
            {
                outer_indexes.insert(iter, i);
                inner_indexes[j]++;
                values.insert(values.begin() + inner, value);
            }
        }
    }

    const T operator()(std::size_t i, std::size_t j) const
    {
        if constexpr (Order == StorageOrder::ROW_WISE)
        {
            auto outer = get_outer_index(i, j);
            auto inner_begin = inner_indexes.begin() + i;
            auto inner_end = inner_indexes.begin() + i + 1;

            auto iter = std::lower_bound(inner_begin, inner_end, outer);

            if (iter != inner_end && *iter == outer)
                return values[std::distance(inner_indexes.begin(), iter)];

            return T{};
        }
        else
        {
            auto inner = get_inner_index(i, j);
            auto outer_begin = outer_indexes.begin() + inner;
            auto outer_end = outer_indexes.begin() + inner + 1;

            auto iter = std::lower_bound(outer_begin, outer_end, i);

            if (iter != outer_end && *iter == i)
                return values[std::distance(outer_indexes.begin(), iter)];

            return T{};
        }
    }

};




#endif // HH_COMPRESSED_MATRIX_HH