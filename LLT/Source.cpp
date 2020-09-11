#include <iostream>
#include "Matrix.cpp"
#define type float
using namespace std;

int main()
{
   Matrix<type> m;
   vector<type> b;
   // Input.
   ifstream in("matrix3.txt");
   m.input(in, b);
   in.close();

   // Solving.
   m.factorization();
   m.direct(b);
   m.reverse(b);

   // Output and clear memory.
   m.~Matrix();
   ofstream out("out.txt");
   for (int i = 0; i < m.getDimension(); i++)
      out << b[i] << " ";
   out.close();
}