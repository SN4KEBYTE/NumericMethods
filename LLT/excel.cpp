#include "excel.h"

void dump_research(ostream &out, const int &k, const vector<float> &xk_float, const vector<double> &xk_double,
    const vector<int> &exact_solution)
{
    out << "k;xk (îäèíàðíàÿ òî÷íîñòü);x* - xk (îäèíàðíàÿ òî÷íîñòü);xk (äâîéíàÿ òî÷íîñòü);x* - xk (äâîéíàÿ òî÷íîñòü);xk (ñêàëÿð. ïðîèçâ.); x* - xk (ñêàëÿð. ïðîèçâ.)\n";

    /*auto s1 = scalar_product(xk_double, xk_double);
    auto s2 = scalar_product(vec_diff(exact_solution, xk_double), vec_diff(exact_solution, xk_double));*/

    for (size_t i = 0; i < xk_float.size(); i++)
    {
        if (i == 0)
            out << k;

        out << ";" << xk_float[i] << ";" << exact_solution[i] - xk_float[i] << ";" << xk_double[i] << ";" <<
            exact_solution[i] - xk_double[i];
    
        //  ÍÅ ÏÎÍÈÌÀÞ, ÊÀÊÈÅ ÑÊÀËßÐÍÛÅ ÏÐÎÈÇÂÅÄÅÍÈß ÒÓÒ ÍÓÆÍÛ
        /*if (i == 0)
            out << ";" << s1 << ";" << s2 << "\n";*/
    }
}