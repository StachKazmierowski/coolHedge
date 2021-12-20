#ifndef HEDGE_CPP_CHOPSTICK_H
#define HEDGE_CPP_CHOPSTICK_H
#include "utils.h"

void test_H_chopstick(int A, int n);

vector<vector<double>> payoff_matrix(int A, int B, int n, int (*aggregation_function)(int, int, int));

vector<int> recurrence_H_chopstick(vector<int> s_A, vector<int> s_B);

vector<int> h_chopstick(vector<int> s_A, vector<int> s_B);

vector<vector<double>> payoff_matrix_chopstick_parallel(int A, int B, int n);

vector<vector<double>> payoff_matrix_chopstick_unsymmetric_parallel(int A, int B, int n);

vector<vector<double>> payoff_matrix_chopstick_symmetric_parallel(int A, int B, int n);


#endif //HEDGE_CPP_CHOPSTICK_H
