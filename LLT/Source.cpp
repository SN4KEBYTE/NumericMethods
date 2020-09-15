#include <iostream>
#include "Matrix.cpp"
#include <string>

using namespace std;

#pragma region Constants, that can be changed

#define type float
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
vector<T>& MatrixDotVector(Matrix<T> m, vector<T> v)
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
   Matrix<type> m;
   vector<type> b;
   // Input.
   ifstream in(pathIn);
   m.input(in, b);
   in.close();

   // Solving.
   m.factorization();
   m.forward(b);
   m.backward(b);

   // Output and clear memory.
   m.~Matrix();
   ofstream out(pathOut);
   ShowVector<type>(b, cout);
   out.close();
}
#pragma endregion