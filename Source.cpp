#include <iostream>
#include "Matrix.cpp"

using namespace std;

template<typename T>
void showMatrix(Matrix<T>);

int main()
{
   ifstream file("matrix.txt", ios_base::in);
   Matrix<float> m(file), L(m.getDimension(), m.getCountElems()), LT(m.getDimension(), m.getCountElems());
   /*showMatrix(m);
   m.s_setElem(4, 1, 11);
   m.s_setElem(4, 1, 6);
   m.s_setElem(3, 3, 6);
   showMatrix(m);*/
   m.s_factorization(&L, &LT);
   showMatrix(&L);
}

template<typename T>
void showMatrix(Matrix<T> *m)
{
   for (size_t i = 0; i < m->getDimension(); i++)
   {
      for (size_t j = 0; j < m->getDimension(); j++)
         cout << m->s_getElem(i, j) << " ";
      cout << endl;
   }
   cout << endl;
}