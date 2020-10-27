#pragma once

#include <vector>
#include <fstream>

const size_t DIAGS_NUM = 9;
const size_t SIDE_DIAGS_NUM = 4;

template<typename T>
class diag_matrix
{
private:
    size_t N = 0, m = 0;
    std::vector<std::vector<T>> diags;
    std::vector<size_t> ig;

    T relative_residual(const std::vector<T> &f, const std::vector<T> &xk)
    {
        return norm(vec_diff(f, dot(xk))) / norm(f);
    }

    T norm(const std::vector<T> x)
    {
        T res = 0;

        for (const auto &el : x)
            res += el * el;

        return sqrt(res);
    }

    std::vector<T> vec_diff(const std::vector<T> &x, const std::vector<T> &y)
    {
        std::vector<T> res(x.size());

        for (size_t i = 0; i < res.size(); i++)
            res[i] = x[i] - y[i];

        return res;
    }

public:
    diag_matrix(std::istream &in)
    {
        diags.resize(DIAGS_NUM);
        ig.resize(SIDE_DIAGS_NUM);

        in >> N >> m;

        for (size_t i = 0; i < diags.size(); i++)
        {
            diags[i].resize(N);

            for (size_t j = 0; j < diags[i].size(); j++)
                in >> diags[i][j];
        }

        ig[0] = 1;

        for (size_t i = 1; i < ig.size(); i++)
            ig[i] = m + i + 1;
    }

    diag_matrix(const size_t &N, const size_t &m, const std::vector<T> diags)
    {
        this->N = N;
        this->m = m;
        this->diags = diags;
    }

    size_t get_dim()
    {
        return N;
    }

    std::vector<T> dot(const std::vector<T> x)
    {
        std::vector<T> y(x.size());

        // processing main diagonal
        for (size_t i = 0; i < N; i++)
        {
            y[i] = diags[4][i] * x[i];

            // lower triangle
            if (i > 0)
                y[i] += diags[3][i] * x[i - ig[0]];

            if (i > ig[1] - 1)
                y[i] += diags[2][i] * x[i - ig[1]];

            if (i > ig[2] - 1)
                y[i] += diags[1][i] * x[i - ig[2]];

            if (i > ig[3] - 1)
                y[i] += diags[0][i] * x[i - ig[3]];

            // upper triangle
            if (i < N - ig[0])
                y[i] += diags[5][i] * x[i + ig[0]];

            if (i < N - ig[1])
                y[i] += diags[6][i] * x[i + ig[1]];

            if (i < N - ig[2])
                y[i] += diags[7][i] * x[i + ig[2]];

            if (i < N - ig[3])
                y[i] += diags[8][i] * x[i + ig[3]];
        }

        return y;
    }

    void step(const double &omega, std::vector<T> &x, std::vector<T> &x_next, std::vector<T> &y,
        const std::vector<T> &f)
    {
        // processing main diagonal separately
        for (size_t i = 0; i < N; i++)
            x_next[i] = diags[4] * x[i];

        // processing upper triangle by rows

        // processing lower triangle by rows

    }

    std::vector<T> jacobi(const double &omega, const std::vector<T> f0, const std::vector<T> &f, const T &eps, const size_t &max_iter)
    {
        std::vector<T> result(N);
    }

    void gauss_seidel()
    {

        return;
    }
};