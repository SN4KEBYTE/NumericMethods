#pragma once

#include <vector>
#include <fstream>

const size_t SIDE_DIAGS_NUM = 4;

template<typename T>
class diag_matrix
{
private:
    size_t N = 0, m = 0;
    std::vector<T> di;
    std::vector<std::vector<T>> ggl(SIDE_DIAGS_NUM);
    std::vector<std::vector<T>> ggu(SIDE_DIAGS_NUM);
    std::vector<size_t> ig(SIDE_DIAGS_NUM);

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
        in >> N >> m;
        di.resize(N);

        for (size_t i = 0; i < di.size(); i++)
            in >> di[i];

        for (size_t i = 0; i < ggl.size(); i++)
        {
            ggl[i].resize(N - 1);

            for (size_t j = 0; j < ggl[i].size(); j++)
                in >> ggl[i][j];
        }

        for (size_t i = 0; i < ggu.size(); i++)
        {
            ggu[i].resize(N - 1);

            for (size_t j = 0; j < ggu[i].size(); j++)
                in >> ggu[i][j];
        }

        ig[0] = 1;
        
        for (size_t i = 1; i < ig.size(); i++)
            ig[i] = m + i;
    }

    diag_matrix(const size_t &N, const size_t &m, const std::vector<T> di, const std::vector<std::vector<T>> &ggl,
        const std::vector<std::vector<T>> &ggu)
    {
        this->N = N;
        this->m = m;
        this->di = di;
        this->ggl = ggl;
        this->ggu = ggu;
    }

    std::vector<T> dot(const std::vector<T> x)
    {
        std::vector<T> y(N);

        for (size_t i = 0; i < di.size(); i++)
            y[i] = di[i] * x[i];

        for (size_t i = 0; i < SIDE_DIAGS_NUM; i++)
            for (size_t j = 0; j < N - ig[i]; j++)
            {
                size_t ir = j + ig[i];
                y[ir] += ggl[i][j] * x[j];
                y[j] += ggu[i][j] * x[ir];
            }

        return y;
    }

    void step(const double &omega, std::vector<T> &x, std::vector<T> &x_next, const std::vector<T> &f)
    {
        // processing main diagonal separately
        for (size_t i = 0; i < di.size(); i++)
            x_next[i] = di[i] * x[i];

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