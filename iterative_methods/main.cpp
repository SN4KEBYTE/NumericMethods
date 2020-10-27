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

#pragma region weight research

const double W0 = 0.05;
const string DATA_PATH = "weight/";

template <typename T>
void dump_table_row(ostream &out, const T &omega, const vector<T> &x, const vector<T> &x_ast, const unsigned &iter)
{
    out << "w;x;x*-x;iterations;\n";
    
    for (size_t i = 0; i < x.size(); i++)
    {
        if (i == 0)
            out << omega;

        out << ";" << x[i] << ";" << x_ast[i] - x[i] << ";";

        if (i == 0)
            out << iter;

        out << ";\n";
    }
}

void weight_research()
{
    ifstream in(DATA_PATH + "matrix.txt");
    diag_matrix<double> a(in);
    in.close();

    in.open(DATA_PATH + "vector0.txt");
    vector<double> f0(a.get_dim());

    for (auto &el : f0)
        in >> el;

    in.close();

    in.open(DATA_PATH + "vector.txt");
    vector<double> f(a.get_dim());

    for (auto &el : f)
        in >> el;

    in.close();

    in.open(DATA_PATH + "exact_solution.txt");
    vector<double> es(a.get_dim());

    for (auto &el : es)
        in >> el;

    in.close();

    in.open(DATA_PATH + "params.txt");
    double eps = 0.;
    unsigned max_iter = 0;

    in >> eps >> max_iter;
    in.close();

    ofstream out("jacobi_res.csv");

    // Jacobi
    for (double w = W0; w <= 1; w += W0)
    {
        unsigned iter = 0;
        auto jacobi_res = a.jacobi(w, f0, f, eps, max_iter, iter);
        dump_table_row(out, w, jacobi_res, es, iter);
    }

    // Gauss-Seidel
    for (double w = W0; w <= 1.95; w += W0)
    {
        unsigned iter = 0;
        auto gs_res = a.gauss_seidel(w, f0, f, eps, max_iter, iter);
        dump_table_row(out, w, gs_res, es, iter);
    }
}

#pragma endregion

int main()
{
    initial_testing();

    system("pause");

    return 0;
}