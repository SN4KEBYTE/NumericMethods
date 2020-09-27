#include <iostream>
#include <limits>
#include "Matrix.cpp"
#include "MatrixGenerator.cpp"
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

template<typename T>
int gauss(vector<vector<T>>& A, vector<T>& x) {
    int n = (int)A.size();
    int m = (int)A[0].size() - 1;

    vector<int> where(m, -1);

    for (int col = 0, row = 0; col < m && row < n; ++col)
    {
        int sel = row;

        for (int i = row; i < n; ++i)
            if (abs(A[i][col]) > abs(A[sel][col]))
                sel = i;

        if (abs(A[sel][col]) < EPS)
            continue;

        for (int i = col; i <= m; ++i)
            swap(A[sel][i], A[row][i]);

        where[col] = row;

        for (int i = 0; i < n; ++i)
            if (i != row)
            {
                double c = A[i][col] / A[row][col];

                for (int j = col; j <= m; ++j)
                    A[i][j] -= A[row][j] * c;
            }
        
        ++row;
    }

    x.assign(m, 0);
    
    for (int i = 0; i < m; ++i)
        if (where[i] != -1)
            x[i] = A[where[i]][m] / A[where[i]][i];

    for (int i = 0; i < n; ++i)
    {
        double sum = 0;

        for (int j = 0; j < m; ++j)
            sum += x[j] * A[i][j];

        if (abs(sum - A[i][m]) > EPS)
            return 0;
    }

    for (int i = 0; i < m; ++i)
        if (where[i] == -1)
            return INFINITY;

    return 1;
}

template<typename T>
void input(ifstream& in, vector<vector<T>>& A)
{
    int dim = 0;
    in >> dim;
    A.resize(dim);

    for (int i = 0; i < dim; i++)
    {
        A[i].resize(dim + 1);

        for (int j = 0; j < dim + 1; j++)
            in >> A[i][j];
    }
}
#pragma endregion

#pragma region Main
int main()
{
    system("pause");
}
#pragma endregion