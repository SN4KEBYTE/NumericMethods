#pragma once

#include <vector>
#include <fstream>
#include <iostream>

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

        for (size_t i = 0; i < N; i++)
        {
            // main diagonal
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

    std::vector<T> step(const double &omega, const std::vector<T> &xk, const std::vector<T> &xk_1, const std::vector<T> &f)
    {
        std::vector<T> y(xk.size());
        T sum1 = 0, sum2 = 0; // sum1 - lower (xk1), sum2 - upper and main diag (xk)

        for (size_t i = 0; i < N; i++)
        {
            // main diagonal
            sum2 += diags[4][i] * xk[i];

            // lower triangle
            if (i > 0)
                sum1 += diags[3][i] * xk_1[i - ig[0]];

            if (i > ig[1] - 1)
                sum1 += diags[2][i] * xk_1[i - ig[1]];

            if (i > ig[2] - 1)
                sum1 += diags[1][i] * xk_1[i - ig[2]];

            if (i > ig[3] - 1)
                sum1 += diags[0][i] * xk_1[i - ig[3]];

            // upper triangle
            if (i < N - ig[0])
                sum2 += diags[5][i] * xk[i + ig[0]];

            if (i < N - ig[1])
                sum2 += diags[6][i] * xk[i + ig[1]];

            if (i < N - ig[2])
                sum2 += diags[7][i] * xk[i + ig[2]];

            if (i < N - ig[3])
                sum2 += diags[8][i] * xk[i + ig[3]];
        }

        for (size_t i = 0; i < N; i++)
            y[i] = xk[i] + omega / diags[4][i] * (f[i] - sum1 - sum2);
        
        return y;
    }

    std::vector<T> jacobi(const double &omega, const std::vector<T> f0, const std::vector<T> &f, const T &eps, const size_t &max_iter)
    {
        auto xk = f0;
        std::vector<T> xk_1;
        unsigned iter = 0;
        T rr = 0;

        do
        {
            xk_1 = step(omega, xk, xk, f);
            rr = relative_residual(f, xk_1);

            std::cout << "Iteration ¹" << iter << "; rr = " << rr << std::endl;

            xk = xk_1;
            iter++;
        } while (rr >= eps && iter < max_iter);

        return xk_1;
    }

    std::vector<T> gauss_seidel(const double &omega, const std::vector<T> f0, const std::vector<T> &f, const T &eps, const size_t &max_iter)
    {
        auto xk = f0;
        std::vector<T> xk_1(xk.size());
        unsigned iter = 0;
        T rr = 0;

        do
        {
            xk_1 = step(omega, xk, xk_1, f);
            rr = relative_residual(f, xk_1);

            std::cout << "Iteration ¹" << iter << "; rr = " << rr << std::endl;

            xk = xk_1;
            iter++;
        } while (rr >= eps && iter < max_iter);

        return xk_1;
    }
};