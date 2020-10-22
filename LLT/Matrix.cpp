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
   std::vector<T>& au = al;
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

   void from_data(const int& dim, const std::vector<T>& al, const std::vector<T>& di, const std::vector<int>& ia)
   {
      this->dim = dim;
      this->al = al;
      this->di = di;
      this->ia = ia;
   }

   void ak_add(T val)
   {
      di[0] += val;
   }
   void display(std::ostream& out)
   {
      for (int i = 0; i < dim; i++)
      {
         for (int j = 0; j < i - (ia[i + 1] - ia[i]); j++)
            out << 0 << " ";
         for (int j = ia[i]; j < ia[i + 1]; j++)
            out << al[j] << " ";
         out << di[i];
         out << std::endl;
      }
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
   void factorization(Matrix<T>& LU)
   {
      for (int i = 0; i < dim; i++)
      {
         // i - index of row
         // j - index of column
         int i0 = ia[i];
         int i1 = ia[i + 1];
         int j = i - (i1 - i0);
         T sum_di = 0;

         for (int k = i0; k < i1; k++, j++)
         {
            int j0 = ia[j];
            int j1 = ia[j + 1];

            T sum = 0;

            int ki = i0;
            int kj = j0;

            int kui = k - i0;
            int kuj = j1 - j0;
            int kur = kui - kuj;

            if (kur > 0)
               ki += kur;
            else
               kj -= kur;

            for (; ki < k; ki++, kj++)
               sum += al[ki] * al[kj];

            LU.al[k] = (al[k] - sum) / LU.di[j];
            sum_di += LU.al[k] * LU.al[k];
         }
         LU.di[i] = sqrt(di[i] - sum_di);
      }
   }

   void factorization_d(Matrix<T>& LU)
   {
      for (int i = 0; i < dim; i++)
      {
         // i - index of row
         // j - index of column
         int i0 = ia[i];
         int i1 = ia[i + 1];
         int j = i - (i1 - i0);
         T sum_di = 0;

         for (int k = i0; k < i1; k++, j++)
         {
            int j0 = ia[j];
            int j1 = ia[j + 1];

            double sum = 0;

            int ki = i0;
            int kj = j0;

            int kui = k - i0;
            int kuj = j1 - j0;
            int kur = kui - kuj;

            if (kur > 0)
               ki += kur;
            else
               kj -= kur;

            for (; ki < k; ki++, kj++)
               sum += al[ki] * al[kj];

            LU.al[k] = (al[k] - sum) / LU.di[j];
            sum_di += LU.al[k] * LU.al[k];
         }
         LU.di[i] = sqrt(di[i] - sum_di);
      }
   }

   void forward(std::vector<T>& y, std::vector<T>& b)
   {
      for (int i = 0; i < dim; i++)
      {
         T elem = b[i];
         int i0 = ia[i];
         int i1 = ia[i + 1];

         // â k îïðåäåëÿåòñÿ ñìåùåíèå
         for (int j = ia[i], k = i1 - i0 < i ? i - (i1 - i0) : 0; j < i1; j++, k++)
            elem -= al[j] * b[k];

         elem /= di[i];
         y[i] = elem;
      }
   }

   void forward_d(std::vector<T>& y, std::vector<T>& b)
   {
      for (int i = 0; i < dim; i++)
      {
         double elem = b[i];
         int i0 = ia[i];
         int i1 = ia[i + 1];

         // â k îïðåäåëÿåòñÿ ñìåùåíèå
         for (int j = ia[i], k = i - (i1 - i0); j < i1; j++, k++)
            elem -= al[j] * b[k];

         elem /= di[i];
         y[i] = elem;
      }
   }

   void backward_modify(std::vector<T>& x, const std::vector<T>& y)
   {
      for (int i = dim - 1; i >= 0; i--)
      {
         int i0 = ia[i];
         int i1 = ia[i + 1];
         int j = i - (i1 - i0);
         T xi = x[i] = y[i] / di[i];

         for (int k = i0; k < i1; k++, j++)
            x[j] -= au[k] * xi;
      }
   }


   void backward(std::vector<T>& x, const std::vector<T>& y)
   {
      for (int i = dim - 1; i >= 0; i--)
      {
         T sum = 0;
         T xi = x[i] = y[i] / di[i];

         for (int j = dim - 1; j > i; j--)
         {
            int bias = ia[j + 1] - ia[j] - j + i;
            if (bias >= 0)
               sum += al[ia[j] + bias] * y[j];
         }
         sum /= di[i];
         x[i] -= sum;
      }
   }

   void backward_d(std::vector<T>& x, const std::vector<T>& y)
   {
      for (int i = dim - 1; i >= 0; i--)
      {
         double sum = 0;
         T xi = x[i] = y[i] / di[i];

         for (int j = dim - 1; j > i; j--)
         {
            int smeshenie = ia[j + 1] - ia[j] - j + i;
            if (smeshenie >= 0)
               sum += al[ia[j] + smeshenie] * y[j];
         }
         sum /= di[i];
         x[i] -= sum;
      }
   }
#pragma endregion

#pragma endregion
};
