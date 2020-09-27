#pragma once

#include <limits>
#include <string>
#include <typeinfo>

#define type double

using namespace std;

// for Gaussian
const double EPS = 1e-14; // FLOAT 1e-6, DOUBLE 1e-14

// correct output precision
constexpr int PREC = numeric_limits<type>::max_digits10;

// get string representation of current type to sort results
const string TYPE_STR = typeid(type).name();