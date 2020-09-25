#include <iostream>
#include "Matrix.cpp"
#include "MatrixGenerator.cpp"
#include <string>

using namespace std;

#pragma region Constants, that can be changed

#define type float
const double EPS = 1e-6; // FLOAT 1e-6, DOUBLE 1e-14
const int numTest = 3;

#pragma endregion

#pragma region Input and output paths

const string pathIn = "test" + to_string(numTest) + "/matrix.txt";
const string pathOut = "test" + to_string(numTest) + "/out.txt";

#pragma endregion

#pragma region Output methods

template<typename T>
void ShowMatrix(Matrix<T>& m, ostream& out)
{
   for (int i = 0; i < m.getDimension(); i++)
   {
      for (int j = 0; j < m.getDimension(); j++)
         out << m.getElem(i, j) << " ";
   
      out << endl;
   }
   out << endl;
}

template<typename T>
void ShowVector(vector<T>& v, ostream& out)
{
   for (size_t i = 0; i < v.size(); i++)
      out << v[i] << " ";
   out << endl;
}

#pragma endregion

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
   for (int col = 0, row = 0; col < m && row < n; ++col) {
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
         if (i != row) {
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
   for (int i = 0; i < n; ++i) {
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
   //ofstream out("matrixAK.txt");
   //m.Gilbert(out, 4);
   //Matrix<type> m;
   //vector<type> b;
   //// Input.
   //ifstream in("matrixAK.txt");//pathIn);
   //m.input(in, b);
   //in.close();

   //// Solving.
   //m.factorization();
   //m.forward(b);
   //m.backward(b);

   //// Output and clear memory.
   //m.~Matrix();
   //ofstream out("out.txt");
   //cout.precision(10);
   //cout.setf(std::ios::fixed);
   //ShowVector<type>(b, cout);
   //out.close();
}
#pragma endregion