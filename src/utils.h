#ifndef HEDGE_CPP_UTILS_H
#define HEDGE_CPP_UTILS_H
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <chrono>
#include <omp.h>

using namespace std;

vector<vector<int>> clash_matrix(vector<int> s_A, vector<int> s_B);

void print_matrix(vector<vector<int>> matrix);

void print_matrix(vector<vector<double>> matrix);

int sign(int a);

bool check_if_contains(vector<vector<int>> vectors, vector<int> vector);

vector<vector<int>> next_divide(vector<vector<int>> divide);

vector<vector<int>> all_permutations(vector<int> s_A);

vector<vector<int>> divides_strategies(int resources, int fields_number);

int k_W(vector<int> s_A, vector<int> s_B);

int k_L(vector<int> s_A, vector<int> s_B);

vector<vector<int>> h(vector<int> s_A, vector<int> s_B);

int factorial(int n);

int matrix_sum(vector<vector<int>> matrix);

vector<int> find_area_vector(vector<vector<int>> clash_mat, int area);

void print_vector(vector<int> v);

vector<vector<int>> empty_matrix(int size);

int height_to_cutoff(vector<int> L, vector<int> T, int current_height);

int width_to_cutoff(vector<int> L, vector<int> T, int current_height);

bool has_single_value(vector<int> v);

int same_values_at_tail(vector<int> v);

vector<int> add_vectors(vector<int> a, vector<int> b);

bool assert_matrices(vector<vector<int>> matrix_1, vector<vector<int>> matrix_2);

bool assert_matrices(vector<vector<double>> matrix_1, vector<vector<double>> matrix_2);

vector<int> trim_vector(vector<int> v, int number_of_elements);

vector<vector<int>> find_knots(vector<int> L, vector<int> T);

int newton_symbol(int n, int k);

int single_type_rectangle(int cols_num, int rows_num, int rooks_num);

int min_number_of_rooks(int i, int j, int n);

vector<vector<int>> recurrence_H(vector<int> s_A, vector<int> s_B);

void test_H(int A, int n);

int blotto(int k_W, int k_L, int n);

int attack(int k_W, int k_L, int n);

int majoritarian(int k_W, int k_L, int n);

int chopstick(int k_W, int k_L, int n);

vector<vector<double>> payoff_matrix(int A, int B, int n, int (*aggregation_function)(int, int, int));

void number_of_different_clash_matrices(int A, int B, int n, int (*aggregation_function)(int, int, int));

double symmetrized_payoff(vector<int> s_A, vector<int> s_B, int (*aggregation_function)(int, int, int));

vector<vector<double>> payoff_matrix_parallel(int A, int B, int n, int (*aggregation_function)(int, int, int));

vector<vector<double>> payoff_matrix_symmetric_parallel(int A, int B, int n, int (*aggregation_function)(int, int, int));

vector<vector<double>> payoff_matrix_unsymmetric_parallel(int A, int B, int n, int (*aggregation_function)(int, int, int));

bool is_symmetric(int (*aggregation_function)(int, int, int), int fields_number);

#endif
