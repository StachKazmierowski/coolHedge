#include "unsymmetrized_matrices.h"
#include "utils.h"
#include <tuple>

vector<vector<int>> all_strategies(int A, int n){
    vector<vector<int>> strategies_symmetrized = divides_strategies(A, n);
    vector<vector<int>> strategies_unsymmetrized = vector<vector<int>>();
    for (vector<int> symm_strategy : strategies_symmetrized){
        for(vector<int> strat : all_permutations(symm_strategy)){
            strategies_unsymmetrized.push_back(strat);
        }
    }
    return strategies_unsymmetrized;
}

vector<vector<int>> unsymmetrized_payoff_matrix_chopstick(int A, int B, int n){
    vector<vector<int>> A_strategies = all_strategies(A, n);
    vector<vector<int>> B_strategies = all_strategies(B, n);
    vector<vector<int>> matrix = vector<vector<int>>();
    for(int i = 0 ; i < A_strategies.size(); i++){
        vector<int> row;
        for(int j = 0 ; j < B_strategies.size(); j++){
            row.push_back(unsymmetrized_payoff_chopstick(A_strategies.at(i), B_strategies.at(j)));
        }
        matrix.push_back(row);
    }
    return matrix;
}

int unsymmetrized_payoff_chopstick(vector<int> strat_A, vector<int> strat_B){
    int result = 0;
    for(int i = 0; i < strat_A.size(); i++){
        result += sign(strat_A.at(i) - strat_B.at(i));
    }
    return sign(result);
}

tuple<vector<vector<double>>, vector<vector<int>>, vector<vector<int>>> payoff_matrix_unsymmetrized_chopstick_parallel(int A, int B, int n) {
    if(A == B){
        cout << "symmetric" << endl;
        return payoff_matrix_unsymmetrized_chopstick_symmetric_parallel(A, B, n);
    }
    cout << "unsymmetric" << endl;
    return payoff_matrix_unsymmetrized_chopstick_unsymmetric_parallel(A, B, n);
}

tuple<vector<vector<double>>, vector<vector<int>>, vector<vector<int>>>  payoff_matrix_unsymmetrized_chopstick_symmetric_parallel(int A, int B, int n){
    vector<vector<double>> matrix;
    vector<vector<int>> A_strategies = all_strategies(A, n);
    for(int i = 0 ; i < A_strategies.size(); i++){
        vector<double> row;
        for(int j = 0 ; j < A_strategies.size(); j++){
            row.push_back(0);
        }
        matrix.push_back(row);
    }
    int i;
    int j;
#pragma omp parallel for default(none) private(i, j) shared(A_strategies, matrix)
    for(i = 0 ; i < A_strategies.size(); i++){
        for(j = i + 1 ; j < A_strategies.size(); j++){
            double payoff =  unsymmetrized_payoff_chopstick(A_strategies.at(i), A_strategies.at(j));
            matrix.at(i).at(j) = payoff;
            matrix.at(j).at(i) = -payoff;
        }
    }
    return tie(matrix, A_strategies, A_strategies);
}

tuple<vector<vector<double>>, vector<vector<int>>, vector<vector<int>>>  payoff_matrix_unsymmetrized_chopstick_unsymmetric_parallel(int A, int B, int n){
    vector<vector<double>> matrix;
    vector<vector<int>> A_strategies = divides_strategies(A, n);
    vector<vector<int>> B_strategies = divides_strategies(B, n);
    for(int i = 0 ; i < A_strategies.size(); i++){
        vector<double> row;
        for(int j = 0 ; j < B_strategies.size(); j++){
            row.push_back(0);
        }
        matrix.push_back(row);
    }
    int i;
    int j;
#pragma omp parallel for default(none) private(i, j) shared(A_strategies, B_strategies, matrix)
    for(i = 0 ; i < A_strategies.size(); i++){
        for(j = 0 ; j < B_strategies.size(); j++){
            matrix.at(i).at(j) = unsymmetrized_payoff_chopstick(A_strategies.at(i), B_strategies.at(j));
        }
    }
    return tie(matrix, A_strategies, B_strategies);
}