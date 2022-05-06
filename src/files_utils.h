#ifndef HEDGE_CPP_FILES_UTILS_H
#define HEDGE_CPP_FILES_UTILS_H
#include "utils.h"
#include "chopstick.h"
#include <fstream>
#include <iomanip>
#include "unsymmetrized_matrices.h"
using namespace std;
static string MATRICES_PATH = "../payoff_matrices_symmetrized/";
static string TIMES_PATH = "../times_symmetrized/";
static string UNSYMMETRIZED_MATRICES_PATH = "../payoff_matrices_unsymmetrized/";
static string UNSYMMETRIZED_TIMES_PATH = "../times_unsymmetrized/";

void save_matrix_to_file(vector<vector<double>> matrix, vector<vector<int>> row_strategies, vector<vector<int>> col_strategies, string path);

void write_col_strategy(ofstream& my_file, vector<int> strategy);

void write_col_strategies(ofstream& my_file, vector<vector<int>> strategies);

void write_rows(ofstream& my_file, vector<vector<double>> matrix, vector<vector<int>> row_strategies);

void save_time_to_file(long int miliseconds, int A, int B, int n, string path);

void find_and_save_chopstick_matrix(int A, int B, int n);

void find_and_save_unsymmetrized_chopstick_matrix(int A, int B, int n);


#endif //HEDGE_CPP_FILES_UTILS_H
