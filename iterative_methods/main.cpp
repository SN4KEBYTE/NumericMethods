#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "diag_matrix.cpp"

using namespace std;

#pragma region initial testing

const string TESTS_PATH = "initial_testing/";
const unsigned INITIAL_TEST_NUM = 1;

void initial_testing()
{
    unsigned dummy_iter = 0;

    ifstream in;
    ofstream out;

    for (unsigned i = 1; i <= INITIAL_TEST_NUM; i++)
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

        auto jacobi_res = a.jacobi(omega, f0, f, eps, max_iter, dummy_iter);
        auto gs_res = a.gauss_seidel(omega, f0, f, eps, max_iter, dummy_iter);

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

const double W0 = 0.01;
const double W_STEP = 0.01;
const double W_JAC_END = 1.;
const double W_GS_END = 1.99;

const string WEIGHT_PATH = "weight/";
const unsigned WEIGHT_TEST_NUM = 2;

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
    for (unsigned i = 1; i <= WEIGHT_TEST_NUM; i++)
    {
        cout << "WEIGHT RESEARCH CASE #" << i << endl;

        string cur_test_path = WEIGHT_PATH + to_string(i) + "/";

        ifstream in(cur_test_path + "matrix.txt");
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

        in.open(cur_test_path + "exact_solution.txt");
        vector<double> es(a.get_dim());

        for (auto &el : es)
            in >> el;

        in.close();

        in.open(cur_test_path + "params.txt");
        double eps = 0.;
        unsigned max_iter = 0;

        in >> eps >> max_iter;
        in.close();

        ofstream out(cur_test_path + "jacobi_res.csv");
        out.imbue(locale(""));

        double w_jac = 0.;
        unsigned iter_jac = max_iter + 1;
        vector<double> jac;

        // Jacobi
        for (unsigned j = 0; j < (W_JAC_END - W0) / W_STEP + 1; j++)
        {
            unsigned iter = 0;
            double w = W0 + j * W_STEP;
            auto jacobi_res = a.jacobi(w, f0, f, eps, max_iter, iter);

            if (iter < iter_jac)
            {
                w_jac = w;
                iter_jac = iter;
                jac = jacobi_res;
            }

            dump_table_row(out, w, jacobi_res, es, iter);
        }

        out.close();
        out.open(cur_test_path + "gauss_seidel_res.csv");

        double w_gs = 0.;
        unsigned iter_gs = max_iter + 1;
        vector<double> gs;

        // Gauss-Seidel
        for (unsigned j = 0; j < (W_GS_END - W0) / W_STEP + 1; j++)
        {
            unsigned iter = 0;
            double w = W0 + j * W_STEP;
            auto gs_res = a.gauss_seidel(w, f0, f, eps, max_iter, iter);

            if (iter < iter_gs)
            {
                w_gs = w;
                iter_gs = iter;
                gs = gs_res;
            }

            dump_table_row(out, w, gs_res, es, iter);
        }

        cout << "W for Jacobi: " << w_jac << endl;
        cout << "Cond for Jacobi: " << a.norm(a.vec_diff(jac, es)) / a.norm(jac) / a.relative_residual(f, jac) << endl;

        cout << "W for Gauss-Seidel: " << w_gs << endl;
        cout << "Cond for Gauss-Seidel: " << a.norm(a.vec_diff(gs, es)) / a.norm(gs) / a.relative_residual(f, gs) << endl;
    }
}

#pragma endregion

int main()
{
    initial_testing();
    weight_research();

    system("pause");

    return 0;
}