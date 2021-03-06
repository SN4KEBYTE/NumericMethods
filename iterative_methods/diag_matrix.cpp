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

    // for solving methods
    std::vector<T> xk;
    std::vector<T> xk_1;

    void step(const double &omega, const std::vector<T> &xk, std::vector<T> &xk_1, std::vector<T> &y, const std::vector<T> &f, 
        T &rr, const bool &is_jacobi = false)
    {
        for (size_t i = 0; i < N; i++)
        {
            // sum1 - lower (xk1), sum2 - upper and main diag (xk)
            // main diagonal
            T sum2 = diags[4][i] * xk[i];
            T sum1 = 0;

            // lower triangle
            if (i > ig[0] - 1)
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

            T tmp = f[i] - sum1 - sum2;
            
            if (is_jacobi)
                rr += tmp * tmp;

            y[i] = xk[i] + omega / diags[4][i] * tmp;
        }

        if (is_jacobi)
            rr = sqrt(rr) / norm(f);
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

        xk.resize(N);
        xk_1.resize(N);
    }

    T relative_residual(const std::vector<T> &f, const std::vector<T> &xk)
    {
        T sum = 0;

        for (size_t i = 0; i < N; i++)
        {
            T elem = diags[4][i] * xk[i];

            if (i > ig[0] - 1)
                elem += diags[3][i] * xk[i - ig[0]];

            if (i > ig[1] - 1)
                elem += diags[2][i] * xk[i - ig[1]];

            if (i > ig[2] - 1)
                elem += diags[1][i] * xk[i - ig[2]];

            if (i > ig[3] - 1)
                elem += diags[0][i] * xk[i - ig[3]];

            if (i < N - ig[0])
                elem += diags[5][i] * xk[i + ig[0]];

            if (i < N - ig[1])
                elem += diags[6][i] * xk[i + ig[1]];

            if (i < N - ig[2])
                elem += diags[7][i] * xk[i + ig[2]];

            if (i < N - ig[3])
                elem += diags[8][i] * xk[i + ig[3]];

            elem = f[i] - elem;
            sum += elem * elem;
        }

        return sqrt(sum) / norm(f);
    }

    T norm(const std::vector<T> &x)
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

    size_t get_dim()
    {
        return N;
    }

    std::vector<T> jacobi(const double &omega, const std::vector<T> &f0, const std::vector<T> &f, const T &eps, const size_t &max_iter,
        unsigned &total_iter)
    {
        xk = f0;
        
        unsigned iter = 0;
        T rr = 0;

        do 
        {
            rr = 0;
            step(omega, xk, xk, xk_1, f, rr, true);

            ++iter;
            std::cout << "Iteration #" << iter << "; rr = " << rr << std::endl;

            xk.swap(xk_1);
        } while (rr >= eps && iter < max_iter);

        total_iter = iter;

        return xk;
    }

    std::vector<T> gauss_seidel(const double &omega, const std::vector<T> &f0, const std::vector<T> &f, const T &eps, const size_t &max_iter,
        unsigned &total_iter)
    {
        xk = f0;
        
        unsigned iter = 0;
        T rr = relative_residual(f, xk);

        while (rr >= eps && iter < max_iter)
        {
            step(omega, xk, xk, xk, f, rr);
            rr = relative_residual(f, xk);

            ++iter;
            std::cout << "Iteration #" << iter << "; rr = " << rr << std::endl;
        }

        total_iter = iter;

        return xk;
    }
};