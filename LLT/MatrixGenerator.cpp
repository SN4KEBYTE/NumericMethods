#include "Matrix.cpp"

template<typename T>
class MatrixGenerator
{
public:
   void Gilbert(std::ofstream out, int dim)
   {
      out << dim;

      // di
      for (int i = 1; i <= dim; i++)
         out << 1 / (2 * i - 1) << " ";
      out << std::endl;

      out << 0 << " ";
      // ia
      for (int i = 0, s = 0; i < dim; i++)
         out << s + i << " ";
      out << std::endl;

      // al
      for (int i = 1; i <= dim; i++)
         for (int j = 1; j < i; j++)
            out << 1 / (i + j - 1);
      out << std::endl;

      // au
      for (int i = 1; i <= dim; i++)
         for (int j = 1; j < i; j++)
            out << 1 / (i + j - 1);
   }
};