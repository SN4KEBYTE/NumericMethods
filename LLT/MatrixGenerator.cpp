#include <vector>
#include <fstream>
#include <random>

class MatrixGenerator
{
public:
   void Gilbert(std::ostream& out, int dim)
   {
      out << dim << std::endl;

      // di
      for (int i = 1; i <= dim; i++)
         out << 1.0 / (2.0 * i - 1) << " ";
      out << std::endl;

      out << 0 << " ";
      // ia
      for (int i = 0, s = 0; i < dim; i++, s += i)
         out << s << " ";
      out << std::endl;

      // al
      for (int i = 1; i <= dim; i++)
         for (int j = 1; j < i; j++)
            out << 1.0 / (i + j - 1.0) << " ";
      out << std::endl;

      // au
      for (int i = 1; i <= dim; i++)
         for (int j = 1; j < i; j++)
            out << 1.0 / (i + j - 1) << " ";
      out << std::endl;
   }

   void Ak(std::ostream& out, int dim, int k)
   {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<> dist(0, 4);

      out << dim << std::endl;
      out.precision(k);
      out.setf(std::ios::fixed);

      std::vector<double> di(dim);
      for (int p = 0; p < dim; p++)
      {
         double elem = 0;
         for (int i = 0; i < dim; i++)
            elem -= dist(gen);
         if (p == 0)
            elem += pow(10, -k);
         out << elem << " ";
         di[p] = elem;
      }
      out << std::endl;

      for (int i = 0; i < dim + 1; i++)
         out << 0 << " ";
      out << std::endl;

      for (int i = 1; i <= dim; i++)
         out << di[i - 1] * i << " ";
      out << std::endl;
   }
};