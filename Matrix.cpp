#include <vector>
#include <fstream>

template<typename T>
class Matrix
{
private:
   int dim = 0;
   std::vector<T> al, au, di;
   std::vector<int> ia, ja;
public:
   int getDimension() { return dim; }

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
      au.resize(size);

      for (int i = 0; i < size; i++)
         inputStream >> al[i];

      for (int i = 0; i < size; i++)
         inputStream >> au[i];

      for (int i = 0; i < dim; i++)
         inputStream >> b[i];
   }

   void factorization()
   {
      di[0] = sqrt(di[0]);
      for (int i = 1; i < dim; i++)
      {
         // l(i, 1)
         T* e = getElemPtr(i, 0);
         if (e != nullptr)
            *e /= di[0];

         // l(i, c)
         // i - индекс текущей строки
         // j - индекс текущего элемента в al (идем по профилю строки)
         for (int j = ia[i] + (e != nullptr ? 1 : 0); j < ia[i + 1]; j++)
         {
            int c = i + j - ia[i + 1]; // столбец в тек. строке
            int wsRow = i - (ia[i + 1] - ia[i]); // количество пустых ячеек в тек. строке
            int wsCol = c - (ia[c + 1] - ia[c]); // количество пустых ячеек в строке колонки

            T r = 0;
            if (wsCol >= wsRow)
            {
               for (int p = ia[c], s = 0; p < ia[c + 1]; p++, s++)
                  r += al[p] * al[ia[i] + s + wsCol - wsRow]; // что-то с индексацией
            }
            else
            {
               for (int p = ia[c] + wsRow - wsCol, s = 0; p < ia[c + 1]; p++, s++)
                  r += al[p] * al[ia[i] + s];
            }

            al[j] -= r;
            al[j] /= di[c];
         }

         // l(i, i)
         for (int p = ia[i]; p < ia[i + 1]; p++)
            di[i] -= al[p] * al[p];
         di[i] = sqrt(di[i]);
      }

      // LT
      au.clear();
      au = *(new std::vector<T>(al));
   }

   void direct(std::vector<T>& b)
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

   void reverse(std::vector<T>& y)
   {
      for (int i = dim - 1; i >= 0; i--)
      {
         T elem = y[i];
         for (int j = ia[dim - i - 1], s = dim - 1; j < ia[dim - i]; j++, s--)
            elem -= au[j] * y[s];
         elem /= di[i];
         y[i] = elem;
      }
   }

   ~Matrix()
   {
      al.clear();
      au.clear();
      di.clear();
      ia.clear();
      ja.clear();
   }
};
