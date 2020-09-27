#pragma once

#include <vector>
#include <fstream>

using namespace std;

template<typename T>
int gauss(vector<vector<T>>& A, vector<T>& x);

template<typename T>
void input(ifstream& in, vector<vector<T>>& A);