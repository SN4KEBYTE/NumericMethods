#pragma once

#include <iostream>
#include <vector>

using namespace std;

void dump_research(ostream &out, const int &k, const vector<float> &xk_float, 
    const vector<float> &xk_float_d, const vector<double> &xk_double,
    const vector<int> &exact_solution);

void dump_gaussian_research(ostream &out, const vector<double> &es, const vector<double> &x_llt,
    const vector<double> &x_gaussian);