#ifndef HEDGE_CPP_FILES_UTILS_H
#define HEDGE_CPP_FILES_UTILS_H
#include "utils.h"
#include <fstream>
#include <iomanip>

using namespace std;

void save_matrix_to_file(vector<vector<double>> matrix, vector<vector<int>> row_strategies, vector<vector<int>> col_strategies, string path);

void write_col_strategy(ofstream& my_file, vector<int> strategy);

void write_col_strategies(ofstream& my_file, vector<vector<int>> strategies);

void write_rows(ofstream& my_file, vector<vector<double>> matrix, vector<vector<int>> row_strategies);


#endif //HEDGE_CPP_FILES_UTILS_H
