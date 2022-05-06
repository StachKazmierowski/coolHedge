#ifndef HEDGE_CPP_UNSYMMETRIZED_MATRICES_H
#define HEDGE_CPP_UNSYMMETRIZED_MATRICES_H
#include <iostream>
#include <vector>
using namespace std;

vector<vector<int>> all_strategies(int A, int n);

vector<vector<int>> unsymmetrized_payoff_matrix_chopstick(int A, int B, int n);

int unsymmetrized_payoff_chopstick(vector<int> strat_A, vector<int> strat_B);

tuple<vector<vector<double>>, vector<vector<int>>, vector<vector<int>>>  payoff_matrix_unsymmetrized_chopstick_symmetric_parallel(int A, int B, int n);

tuple<vector<vector<double>>, vector<vector<int>>, vector<vector<int>>>  payoff_matrix_unsymmetrized_chopstick_unsymmetric_parallel(int A, int B, int n);

#endif //HEDGE_CPP_UNSYMMETRIZED_MATRICES_H
