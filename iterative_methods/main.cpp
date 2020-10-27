#include <iostream>
#include <fstream>
#include <vector>
#include "diag_matrix.cpp"

using namespace std;

int main()
{
    ifstream in("initial_testing/1/matrix.txt");

    diag_matrix<double> m(in);
    vector<double> x(m.get_dim());

    for (size_t i = 0; i < x.size(); i++)
        x[i] = i + 1;

    auto dot_res = m.dot(x);

    for (const auto &el : dot_res)
        cout << el << " ";

    system("pause");

    return 0;
}