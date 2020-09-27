#include "out.h"

template<typename T>
void print_matrix(Matrix<T>& m, ostream& out)
{
    for (int i = 0; i < m.getDimension(); i++)
    {
        for (int j = 0; j < m.getDimension(); j++)
            out << m.getElem(i, j) << " ";

        out << endl;
    }
}

template<typename T>
void print_vector(vector<T>& v, ostream& out)
{
    for (size_t i = 0; i < v.size(); i++)
        out << v[i] << endl;
}