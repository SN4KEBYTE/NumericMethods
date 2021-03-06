#include "excel.h"


void dump_research(ostream &out, const int &k, const vector<float> &xk_float, const vector<float> &xk_float_d, const vector<double> &xk_double,
    const vector<int> &exact_solution)
{
    out << "k;xk (��������� ��������);x* - xk (��������� ��������);xk (������� ��������);x* - xk (������� ��������);xk (������. ������.); x* - xk (������. ������.);\n";

    for (size_t i = 0; i < xk_float.size(); i++)
    {
        if (i == 0)
            out << k;

        out << ";" << xk_float[i] << ";" << exact_solution[i] - xk_float[i] << ";" << xk_double[i] << ";" <<
            exact_solution[i] - xk_double[i] << ";" << xk_float_d[i] << ";" << exact_solution[i] - xk_float_d[i] << ";\n";
    }
}

void dump_gaussian_research(ostream &out, const vector<double> &es, const vector<double> &x_llt,
    const vector<double> &x_gaussian)
{
    out << "x*;x_llt;x_gaussian;x* - x_llt;x* - x_gaussian;\n";

    for (size_t i = 0; i < es.size(); i++)
        out << es[i] << ";" << x_llt[i] << ";" << x_gaussian[i] << ";" <<
            es[i] - x_llt[i] << ";" << es[i] - x_gaussian[i] << ";\n";
}