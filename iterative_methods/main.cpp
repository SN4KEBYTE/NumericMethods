#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "diag_matrix.cpp"

using namespace std;

#pragma region initial testing

const string TESTS_PATH = "initial_testing/";
const unsigned TEST_NUM = 1;

void initial_testing()
{
    ifstream in;
    ofstream out;

    for (unsigned i = 1; i <= TEST_NUM; i++)
    {
        auto cur_test_path = TESTS_PATH + to_string(i) + "/";
        
        in.open(cur_test_path + "matrix.txt");
        diag_matrix<double> a(in);
        in.close();

        in.open(cur_test_path + "vector0.txt");
        vector<double> f0(a.get_dim());

        for (auto &el : f0)
            in >> el;

        in.close();

        in.open(cur_test_path + "vector.txt");
        vector<double> f(a.get_dim());

        for (auto &el : f)
            in >> el;

        in.close();

        in.open(cur_test_path + "params.txt");
        double omega = 0., eps = 0.;
        unsigned max_iter = 0;

        in >> omega >> eps >> max_iter;
        in.close();

        auto jacobi_res = a.jacobi(omega, f0, f, eps, max_iter);
        auto gs_res = a.gauss_seidel(omega, f0, f, eps, max_iter);

        out.open(cur_test_path + "res_jacobi.txt");

        for (const auto &el : jacobi_res)
            out << el << " ";

        out.close();

        out.open(cur_test_path + "res_gauss_seidel.txt");

        for (const auto &el : gs_res)
            out << el << " ";

        out.close();
    }
}

#pragma endregion


int main()
{
    initial_testing();

    system("pause");

    return 0;
}