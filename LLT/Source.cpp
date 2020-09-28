#include <iostream>
#include <limits>
#include "Matrix.cpp"
#include "MatrixGenerator.cpp"
#include "global_const.h"
#include <string>
#include <typeinfo>

using namespace std;

#pragma region Operations with matrix

template<typename T>
vector<T>& MatrixDotVector(Matrix<T>& m, vector<T>& v)
{
    vector<T> res(m.getDimension());

    for (int i = 0; i < m.getDimension(); i++)
        for (int j = 0; j < m.getDimension(); j++)
            res[i] += m.getElem(i, j) * v[j];

    return *(new vector<T>(res));
}
#pragma endregion

#pragma region Main
int main()
{
    system("pause");
}
#pragma endregion