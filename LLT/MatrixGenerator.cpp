#pragma once

#include <vector>
#include <fstream>
#include <random>
#include "Matrix.cpp"

template<typename T, typename V>
std::vector<T> MatrixDotVector(Matrix<T> &m, std::vector<V> &v)
{
    std::vector<T> res(m.getDimension());

    for (int i = 0; i < m.getDimension(); i++)
        for (int j = 0; j < m.getDimension(); j++)
            res[i] += m.getElem(i, j) * v[j];

    return res;
}

class MatrixGenerator
{
public:
    void Gilbert(std::ostream &out, const int &dim)
    {
        Matrix<double> m;
        std::vector<double> diag(dim), al(dim * dim);
        std::vector<int> ia(dim + 1);

        out << dim << std::endl;

        // di
        for (int i = 1; i <= dim; i++)
        {
            double v = 1.0 / (2.0 * i - 1);

            diag[i - 1] = v;
            out << v << " ";
        }

        out << std::endl;

        ia[0] = 0;
        out << 0 << " ";
        // ia
        for (int i = 0, s = 0; i < dim; i++, s += i)
        {
            ia[i + 1] = s;
            out << s << " ";
        }
        out << std::endl;

        // al
        for (int i = 1; i <= dim; i++)
            for (int j = 1; j < i; j++)
            {
                double v = 1.0 / (i + j - 1.0);

                al[(i - 1) * dim + j - 1] = v;
                out << v << " ";
            }

        al.shrink_to_fit();

        out << std::endl;

        m.from_data(dim, al, diag, ia);

        std::vector<double> es(dim);

        for (size_t i = 0; i < dim; i++)
            es[i] = i + 1;

        for (const auto &el : MatrixDotVector(m, es))
            out << el << " ";
    }
};