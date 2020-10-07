#pragma once

#include <vector>
#include <fstream>

template<typename T>
class Matrix
{
#pragma region Private
private:
    int dim = 0;
    std::vector<T> al, di;
    std::vector<T> &au = al;
    std::vector<int> ia;
#pragma endregion

#pragma region Public
public:
    ~Matrix()
    {
        al.clear();
        di.clear();
        ia.clear();
    }

#pragma region Methods
    inline int getDimension() { return dim; }

    T getElem(int i, int j)
    {
        if (i == j)
            return di[i];
        if (i > j)
            if (j < i - ia[i + 1] + ia[i])
                return 0;
            else
                return al[j - i + ia[i + 1]];
        else
            if (i < j - ia[j + 1] + ia[j])
                return 0;
            else
                return au[i - j + ia[j + 1]];
    }

    T* getElemPtr(int i, int j)
    {
        if (i == j)
            return &di[i];
        if (i > j)
            if (j < i - ia[i + 1] + ia[i])
                return nullptr;
            else
                return &al[j - i + ia[i + 1]];
        else
            if (i < j - ia[j + 1] + ia[j])
                return nullptr;
            else
                return &au[i - j + ia[j + 1]];
    }

    void input(std::ifstream& inputStream, std::vector<T>& b)
    {
        inputStream >> dim;
        di.resize(dim);
        ia.resize(dim + 1);
        b.resize(dim);

        for (int i = 0; i < dim; i++)
            inputStream >> di[i];

        for (int i = 0; i < dim + 1; i++)
            inputStream >> ia[i];

        auto size = ia.back();
        al.resize(size);

        for (int i = 0; i < size; i++)
            inputStream >> al[i];

        for (int i = 0; i < dim; i++)
            inputStream >> b[i];
    }
#pragma endregion

#pragma region Solving methods
    // ((dim - 1)^3 * dim^4 / 8)
    void factorization()
    {
        di[0] = sqrt(di[0]);
        
        for (int i = 1; i < dim; i++)
        {
#pragma region l(i, 1)
            T* e = getElemPtr(i, 0);
            if (e != nullptr)
                *e /= di[0];
#pragma endregion

#pragma region l(i, c)
            // i - index of row
            // j - index this element in al (go along the profile)
            for (int j = ia[i] + (e != nullptr ? 1 : 0); j < ia[i + 1]; j++)
            {
                int c = i + j - ia[i + 1];           // column in row
                int wsRow = i - (ia[i + 1] - ia[i]); // the number of whitespaces in the current row
                int wsCol = c - (ia[c + 1] - ia[c]); // the number of whitespaces in row of the column

                T r = 0;
                
                if (wsCol >= wsRow)
                    for (int p = ia[c], s = 0; p < ia[c + 1]; p++, s++)
                        r += al[p] * al[ia[i] + s + wsCol - wsRow];
                else
                    for (int p = ia[c] + wsRow - wsCol, s = 0; p < ia[c + 1]; p++, s++)
                        r += al[p] * al[ia[i] + s];

                al[j] -= r;
                al[j] /= di[c];
            }
#pragma endregion

#pragma region l(i, i)
            for (int p = ia[i]; p < ia[i + 1]; p++)
                di[i] -= al[p] * al[p];
            
            di[i] = sqrt(di[i]);
#pragma endregion
        }
    }

void factorization_d()
    {
        di[0] = sqrt(di[0]);
        
        for (int i = 1; i < dim; i++)
        {
#pragma region l(i, 1)
            T* e = getElemPtr(i, 0);
            if (e != nullptr)
                *e /= di[0];
#pragma endregion

#pragma region l(i, c)
            // i - index of row
            // j - index this element in al (go along the profile)
            for (int j = ia[i] + (e != nullptr ? 1 : 0); j < ia[i + 1]; j++)
            {
                int c = i + j - ia[i + 1];           // column in row
                int wsRow = i - (ia[i + 1] - ia[i]); // the number of whitespaces in the current row
                int wsCol = c - (ia[c + 1] - ia[c]); // the number of whitespaces in row of the column

                double r = 0;
                
                if (wsCol >= wsRow)
                    for (int p = ia[c], s = 0; p < ia[c + 1]; p++, s++)
                        r += al[p] * al[ia[i] + s + wsCol - wsRow];
                else
                    for (int p = ia[c] + wsRow - wsCol, s = 0; p < ia[c + 1]; p++, s++)
                        r += al[p] * al[ia[i] + s];

                al[j] -= r;
                al[j] /= di[c];
            }
#pragma endregion

#pragma region l(i, i)
            double r = 0;
            for (int p = ia[i]; p < ia[i + 1]; p++)
                r -= al[p] * al[p];
            
            di[i] -= r;
            di[i] = sqrt(di[i]);
#pragma endregion
        }
    }


    // (dim^2(dim - 1) / 2)
    void forward(std::vector<T>& b)
    {
        for (int i = 0; i < dim; i++)
        {
            T elem = b[i];
            int d = ia[i + 1] - ia[i];

            for (int j = ia[i], s = d < i ? i - d : 0; j < ia[i + 1]; j++, s++)
                elem -= al[j] * b[s];

            elem /= di[i];
            b[i] = elem;
        }
    }

    // (dim^2(dim - 1) / 2)
    void backward(std::vector<T>& y)
    {
        y[dim - 1] /= di[dim - 1];

        for (int i = dim - 1; i > 0; i--)
        {
            for (int j = ia[i + 1] - 1, s = i - 1; j >= ia[i]; j--, s--)
                y[s] -= au[j] * y[i];
            
            y[i - 1] /= di[i - 1];
        }
    }
#pragma endregion

#pragma endregion
};
