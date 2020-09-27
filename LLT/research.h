#pragma once

#include <fstream>
#include <string>
#include <vector>
#include "global_const.h"
#include "Matrix.cpp"
#include "out.h"
#include "excel.h"

using namespace std;

#pragma region amount of tests and researches
// for initial tests
const int TESTS_NUM = 5;

// for Ak matrices research
const int AK_DIM = 10;
const int AK_NUM = 15;

// for Gilbert matrices research
const int GILBERT_NUM = 9;

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

#pragma endregion

void initial_testing();

void ak_research();

void gilbert_research();

vector<int> get_exact_solution(const int &N);