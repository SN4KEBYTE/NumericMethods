#include <vector>
#include <fstream>
#include <tuple>

using namespace std;

template<typename T>
class Matrix
{
private:
   size_t countElements;
   size_t dimension;

   // Diagonal.
   std::vector<T> di;

   // Lover triangle.
   std::vector<T> al;

   // Upper triangle.
   std::vector<T> au;

   // Indices of rows in array al.
   std::vector<uint32_t> ia;

   // Indices of rows in array au.
   std::vector<uint32_t> ja;

   T* s_getPtrElem(uint32_t row, uint32_t col)
   {
      if (row == col)
         return &di[row];

      if (col < row - ia[row + 1] + ia[row]
         || row < col - ia[col + 1] + ia[col])
         return 0;

      return &al[row > col ? ia[row + 1] - row + col : ia[col + 1] - col + row];
   }

public:
   // Symmetric profile matrix.
   Matrix(size_t dim, size_t countElems)
   {
      dimension = dim;
      countElements = countElems;
      di.resize(dim);
      al.resize(countElems);
      ia.resize(dim + 1);
   }

   // Symmetric profile matrix.
   Matrix(ifstream& input)
   {
      input >> dimension >> countElements;
      di.resize(dimension);
      ia.resize(dimension + 1);
      al.resize(countElements);

      for (size_t i = 0; i < dimension + 1; i++)
         input >> ia[i];
      for (size_t i = 0; i < dimension; i++)
         input >> di[i];
      for (size_t i = 0; i < countElements; i++)
         input >> al[i];
   }

   ~Matrix()
   {
      di.clear();
      al.clear();
      au.clear();
      ia.clear();
      ja.clear();
   }

   size_t getDimension()
   {
      return dimension;
   }

   size_t getCountElems()
   {
      return countElements;
   }

   // For symmetric matrix.
   T s_getElem(uint32_t row, uint32_t col)
   {
      if (row == col)
         return di[row];

      if (col < row - ia[row + 1] + ia[row]
         || row < col - ia[col + 1] + ia[col])
         return 0;

      return al[row > col ? ia[row + 1] - row + col : ia[col + 1] - col + row];
   }

   void s_setElem(uint32_t row, uint32_t col, T elem)
   {
      if (row == col)
         di[row] = elem;
      else
      {
         T* e = s_getPtrElem(row, col);
         if (e == 0)
         {
            al.resize(++countElements);
            for (size_t i = ia[row + 1]; i > ia[row]; i--)
               al[i] = al[i - 1];
            al[ia[row]] = elem;
            ia[row + 1]++;
         }
         else
            *e = elem;
      }
   }

   void s_factorization(Matrix<T>* L, Matrix<T>* LT)
   {
      L->s_setElem(1, 1, sqrt(s_getElem(1, 1)));
      LT->s_setElem(1, 1, sqrt(s_getElem(1, 1)));
      for (size_t j = 1; j < dimension; j++)
      {
         T elem = s_getElem(j, 1) / L->s_getElem(1, 1);
         L->s_setElem(j, 1, elem);
         LT->s_setElem(1, j, elem);
      }

      for (size_t i = 1; i < dimension; i++)
      {
         T elem = s_getElem(i, i);
         for (size_t j = 0; j < i - 1; j++)
         {
            T p = L->s_getElem(i, j);
            elem -= p * p;
         }
         elem = sqrt(elem);
         L->s_setElem(i, i, elem);
         LT->s_setElem(i, i, elem);
      }

      for (size_t i = 1; i < dimension - 1; i++)
      {
         for (size_t j = i + 1; j < dimension; j++)
         {
            T elem = s_getElem(j, i);
            for (size_t p = 0; p < i - 1; p++)
               elem -= L->s_getElem(i, p) * L->s_getElem(j, p);
            elem /= L->s_getElem(i, i);
            L->s_setElem(j, i, elem);
            LT->s_setElem(i, j, elem);
         }
      }
   }
};
