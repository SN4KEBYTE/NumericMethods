#pragma once

#include <vector>

using namespace std;

template<typename T, typename S>
T scalar_product(const vector<T> &a, const vector<S> &b)
{
    T res = 0;

    for (size_t i = 0; i < a.size(); i++)
        res += a[i] * (T)b[i];

    return res;
}

template<typename T, typename S>
vector<T> vec_diff(const vector<T> &a, const vector<S> &b)
{
    vector<T> res(a.size());

    for (size_t i = 0; i < res.size(); i++)
        res[i] = a[i] - (T)b[i];

    return res;
}
