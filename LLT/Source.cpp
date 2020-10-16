#include <iostream>
#include <limits>
#include "MatrixGenerator.cpp"
#include "global_const.h"
//#include "research.h"

using namespace std;

int main()
{
   Matrix<double> m;
   ifstream in("ak/tests/ak.txt");
   vector<double> b;
   m.input(in, b);
   m.factorization(m);
   m.display(cout);
   m.forward(b, b);
   m.backward(b, b);
   for (int i = 0; i < 10; i++)
      cout << b[i] << ' ';
    //// --- GILBERT ---
    //MatrixGenerator gen;
    //
    //for (int i = 2; i <= GILBERT_NUM; i++)
    //{
    //    ofstream out(GILBERT_TESTS_PATH + "gilbert" + to_string(i) + ".txt");
    //    gen.Gilbert(out, i);
    //}

    //gilbert_research();

    //// --- AK ---
    //ak_research();

    //// --- GAUSSIAN ---
    ////gaussian_research();
}