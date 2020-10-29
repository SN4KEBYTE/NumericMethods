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

const double W0 = 0.1;
const double W_STEP = 0.1;
const double W_JAC_END = 1.;
const double W_GS_END = 1.99;

const string WEIGHT_PATH = "weight/";
const unsigned WEIGHT_TEST_NUM = 2;

template <typename T>
void dump_table_row(ostream &out, const T &omega, const vector<T> &x, const vector<T> &x_ast, const unsigned &iter, const T &cond)
{
    out << "w;x;x*-x;iterations;cond;\n";
    
    for (size_t i = 0; i < x.size(); i++)
    {
        if (i == 0)
            out << omega;

        out << ";" << x[i] << ";" << x_ast[i] - x[i] << ";";

        if (i == 0)
            out << iter << ";" << cond;

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

#pragma region Jacobi Global Optimal Weight
        ofstream out(cur_test_path + "jacobi_res_global.csv");
        out.imbue(locale(""));

        double w_jac = 0.;
        unsigned iter_jac = max_iter + 1;
        vector<double> jac;
        double cond_jac = 0;

        // Jacobi
        // search for global optimal weight
        for (unsigned j = 0; j < (W_JAC_END - W0) / W_STEP + 1; j++)
        {
            unsigned iter = 0;
            double w = W0 + j * W_STEP;
            auto jacobi_res = a.jacobi(w, f0, f, eps, max_iter, iter);
            double cond = a.norm(a.vec_diff(jacobi_res, es)) / a.norm(jacobi_res) / a.relative_residual(f, jacobi_res);

            if (iter < iter_jac)
            {
                w_jac = w;
                iter_jac = iter;
                jac = jacobi_res;
                cond_jac = cond;
            }

            dump_table_row(out, w, jacobi_res, es, iter, cond);
        }

        cout << "GLOBAL" << endl;
        cout << "W for Jacobi: " << w_jac << endl;
        cout << "Cond for Jacobi: " << cond_jac << endl;
        cout << "Iter for Jacobi: " << iter_jac << endl;

        out.close();
#pragma endregion        
     
#pragma region Jacobi Local Optimal Weight
        // Jacobi
        // search for local optimal weight
        out.open(cur_test_path + "jacobi_res_local.csv");
        
        double local_w_jac = 0;
        unsigned local_iter_jac = max_iter + 1;
        vector<double> local_jac;
        double local_cond_jac = 0;

        for (unsigned j = 1; j <= 19; j++)
        {
            unsigned iter = 0;
            double w = w_jac + 0.01 * ((int)j - 10);

            if (w > 1.)
                continue;

            auto jacobi_res = a.jacobi(w, f0, f, eps, max_iter, iter);
            double cond = a.norm(a.vec_diff(jacobi_res, es)) / a.norm(jacobi_res) / a.relative_residual(f, jacobi_res);

            if (iter < local_iter_jac)
            {
                local_w_jac = w;
                local_iter_jac = iter;
                local_jac = jacobi_res;
                local_cond_jac = cond;
            }

            dump_table_row(out, w, jacobi_res, es, iter, cond);
        }

        cout << "LOCAL" << endl;
        cout << "W for Jacobi: " << local_w_jac << endl;
        cout << "Cond for Jacobi: " << local_cond_jac << endl;
        cout << "Iter for Jacobi: " << iter_jac << endl;

        out.close();
#pragma endregion

#pragma region Gauss-Seidel Global Optimal Weight
        out.open(cur_test_path + "gauss_seidel_res_global.csv");

        double w_gs = 0.;
        unsigned iter_gs = max_iter + 1;
        vector<double> gs;
        double cond_gs = 0;

        // Gauss-Seidel
        for (unsigned j = 0; j < (W_GS_END - W0) / W_STEP + 1; j++)
        {
            unsigned iter = 0;
            double w = W0 + j * W_STEP;
            auto gs_res = a.gauss_seidel(w, f0, f, eps, max_iter, iter);
            double cond = a.norm(a.vec_diff(gs_res, es)) / a.norm(gs_res) / a.relative_residual(f, gs_res);

            if (iter < iter_gs)
            {
                w_gs = w;
                iter_gs = iter;
                gs = gs_res;
                cond_gs = cond;
            }

            dump_table_row(out, w, gs_res, es, iter, cond);
        }

        cout << "GLOBAL" << endl;
        cout << "W for Gauss-Seidel: " << w_gs << endl;
        cout << "Cond for Gauss-Seidel: " << cond_gs << endl;
        cout << "Iter for Gauss-Seidel: " << iter_gs << endl;

        out.close();
#pragma endregion

#pragma region Gauss-Seidel Local Optimal Weight
        out.open(cur_test_path + "gauss_seidel_res_local.csv");

        double local_w_gs = 0.;
        unsigned local_iter_gs = max_iter + 1;
        vector<double> local_gs;
        double local_cond_gs = 0;

        // Gauss-Seidel
        for (unsigned j = 1; j <= 19; j++)
        {
            unsigned iter = 0;
            double w = w_gs + 0.01 * ((int)j - 10);

            if (w >= 2.)
                continue;

            auto gs_res = a.gauss_seidel(w, f0, f, eps, max_iter, iter);
            double cond = a.norm(a.vec_diff(gs_res, es)) / a.norm(gs_res) / a.relative_residual(f, gs_res);

            if (iter < local_iter_gs)
            {
                local_w_gs = w;
                local_iter_gs = iter;
                local_gs = gs_res;
                local_cond_gs = cond;
            }

            dump_table_row(out, w, gs_res, es, iter, cond);
        }

        cout << "LOCAL" << endl;
        cout << "W for Gauss-Seidel: " << local_w_gs << endl;
        cout << "Cond for Gauss-Seidel: " << local_cond_gs << endl;
        cout << "Iter for Gauss-Seidel: " << local_iter_gs << endl;

        out.close();
#pragma endregion
    }
}

#pragma endregion

int main()
{
    weight_research();

    system("pause");

    return 0;
}