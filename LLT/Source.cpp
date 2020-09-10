#include <iostream>
#include "Matrix.cpp"
#define type float
using namespace std;

int main()
{
   Matrix<type> m;
   vector<type> b;
   ifstream in("matrix3.txt");
   m.input(in, b);
   m.factorization();
   m.direct(b);
   m.reverse(b);
   m.~Matrix();
   for (int i = 0; i < m.getDimension(); i++)
      cout << b[i] << " ";
}