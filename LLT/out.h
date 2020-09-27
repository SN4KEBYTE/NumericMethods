#pragma once

#include <vector>
#include <iostream>
#include "Matrix.cpp"

using namespace std;

template<typename T>
void print_matrix(Matrix<T> &m, ostream &out);

template<typename T>
void print_vector(vector<T> &v, ostream &out);