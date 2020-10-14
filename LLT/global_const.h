#pragma once

#include <limits>
#include <string>
#include <typeinfo>

#define type double

using namespace std;

// correct output precision
constexpr int PREC = numeric_limits<type>::max_digits10;

// get string representation of current type to sort results
const string TYPE_STR = typeid(type).name();