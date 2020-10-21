#pragma once

#include <fstream>
#include <string>
#include <vector>
#include "Matrix.cpp"
#include "excel.h"

using namespace std;

#pragma region amount of tests and researches
// for initial tests
const int TESTS_NUM = 5;

// for Ak matrices research
const int AK_DIM = 10;
const int AK_NUM = 15;

// for Gilbert matrices research
const int GILBERT_NUM = 15;

// for Gaussian research
const int G_DIM = 10;
const double EPS = 1e-14; // FLOAT 1e-6, DOUBLE 1e-14

#pragma endregion

#pragma region input and output paths
// for initial tests
const string INITIAL_TESTS_PATH = "initial/tests/";
const string INITIAL_RESULTS_PATH = "initial/results/";

// for Ak matrices research
const string AK_TESTS_PATH = "ak/tests/";
const string AK_RESULTS_PATH = "ak/results/";

// for Gilbert matrices research
const string GILBERT_TESTS_PATH = "gilbert/tests/";
const string GILBERT_RESULTS_PATH = "gilbert/results/";

// for Gaussian research
const string GAUSSIAN_TESTS_PATH = "gaussian/tests/";
const string GAUSSIAN_RESULTS_PATH = "gaussian/results/";

#pragma endregion

void initial_testing();

void ak_research();

void gilbert_research();

void gaussian_research();

template <typename T>
int gauss(std::vector<std::vector<T>> &A, std::vector<T> &x);

template<typename T>
vector<T> get_exact_solution(const int &N);

template<typename T>
void print_vector(std::vector<T> &v, std::ostream &out);