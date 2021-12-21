#include "files_utils.h"

void save_matrix_to_file(vector<vector<double>> matrix, vector<vector<int>> row_strategies, vector<vector<int>> col_strategies, string path){
    ofstream my_file;
    my_file.open (path);
    write_col_strategies(my_file, col_strategies);
    write_rows(my_file, matrix, row_strategies);
    my_file.close();
}

void write_col_strategy(ofstream& my_file, vector<int> strategy){
    my_file << ",[";
    for (int res : strategy){
        my_file << res << ". ";
    }
    my_file << "]";
}

void write_row_strategy(ofstream& my_file, vector<int> strategy){
    my_file << "[";
    for (int res : strategy){
        my_file << res << ". ";
    }
    my_file << "]";
}

void write_col_strategies(ofstream& my_file, vector<vector<int>> strategies){
    for(vector<int> strategy : strategies)
        write_col_strategy(my_file, strategy);
    my_file << "\n";
}

void write_row(ofstream& my_file, vector<int> strategy, vector<double> row){
    write_row_strategy(my_file, strategy);
    for(double value : row)
        my_file << "," << value;
    my_file << "\n";
}

void write_rows(ofstream& my_file, vector<vector<double>> matrix, vector<vector<int>> row_strategies){
    my_file << setprecision(15);
    for(int i = 0; i < row_strategies.size(); i++)
        write_row(my_file, row_strategies[i], matrix[i]);
}
