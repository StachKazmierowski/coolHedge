#include "chopstick.h"

vector<vector<vector<long int>>> empty_H_chopstick(int knots_number, int fields_number){
    vector<vector<vector<long int>>> values; // knot, rooks_num, value
    for(int i = 0; i < knots_number; i++){
        vector<vector<long int>> values_0;
        for(int j = 0; j < fields_number + 1; j++){
            vector<long int> values_1;
            for(int j = 0; j < 2 * fields_number + 1; j++){
                values_1.push_back(0);
            }
            values_0.push_back(values_1);
        }
        values.push_back(values_0);
    }
    return values;
}

int vector_sum(vector<long int> v){
    int result = 0;
    for(int i = 0 ; i < v.size(); i++)
        result += v.at(i);
    return result;
}

vector<long int> h_chopstick(vector<int> s_A, vector<int> s_B){
    int fields_number = s_A.size();
    vector<long int> values;
    for (int i  = 0 ; i <= 2 * fields_number; i++)
        values.push_back(0);
    vector<vector<int>> permutations;
    permutations = all_permutations(s_A);
    for (int i = 0; i < permutations.size(); i++){
        vector<int> permutation = permutations.at(i);
        long int kW = k_W(permutation, s_B);
        long int kL = k_L(permutation, s_B);
        values.at(fields_number + kW - kL)++;
    }
    long int scale = (long int)(factorial(fields_number) / vector_sum(values));
    for(int i = 0; i < values.size(); i++){
        values.at(i) *= scale;
    }
    return values;
}

vector<long int> recurrence_H_chopstick(vector<int> s_A, vector<int> s_B){
    int fields_number = s_A.size();
    vector<vector<int>> clash_mat = clash_matrix(s_A, s_B);
    vector<int> L = find_area_vector(clash_mat, -1);
    vector<int> T = find_area_vector(clash_mat, 0);
    vector<vector<int>> knots = find_knots(L, T);
    vector<vector<vector<long int>>> values = empty_H_chopstick(knots.size(), fields_number);
    vector<int> knot = knots.at(0);
    int i = knot.at(0);
    int j = knot.at(1);
    if(L.at(0) > 0){
        for (int num_rooks = 0; num_rooks <= min(i,j); num_rooks++)
            values[0][num_rooks][fields_number - num_rooks] = single_type_rectangle(i, j, num_rooks);
    } else if (T.at(0) > 0){
        for (int num_rooks = 0; num_rooks <= min(i,j); num_rooks++)
            values[0][num_rooks][fields_number] = single_type_rectangle(i, j, num_rooks);
    } else {
        for (int num_rooks = 0; num_rooks <= min(i,j); num_rooks++)
            values[0][num_rooks][fields_number + num_rooks] = single_type_rectangle(i, j, num_rooks);
    }
    for(int knot_index = 1; knot_index < knots.size() - 1; knot_index++){
        knot = knots.at(knot_index);
        i = knot.at(0);
        j = knot.at(1);
        vector<int> previous_knot = knots[knot_index-1];
        int old_i = previous_knot[0];
        int old_j = previous_knot[1];
        int delta_i = i - old_i;
        int delta_j = j - old_j;
        if(knot_index == 1 && L[0] > 0 && T[0] > 0){ // L and T stripes
            int maximum_rooks_in_T = min(delta_i, j);
            for(int num_rooks = 0; num_rooks <= min(i,j); num_rooks++){
                for(int r_T = 0; r_T <= min(maximum_rooks_in_T, num_rooks); r_T++){
                    int rooks_left = num_rooks - r_T;
                    long int H_tmp = values[0][rooks_left][fields_number - rooks_left];
                    long int bottom = single_type_rectangle(delta_i, j - rooks_left, r_T);
                    values[knot_index][num_rooks][fields_number - rooks_left] = H_tmp * bottom;
                }
            }
        } else if (knot_index == 1 && T[0] == 0 && T[j-1] > 0 && L[0] == 0 && has_single_value(trim_vector(L, j))){ //W and T stripes
            int maximum_rooks_in_T = min(i, delta_j);
            for(int num_rooks = 0; num_rooks <= min(i,j); num_rooks++){
                for(int r_T = 0; r_T <= min(maximum_rooks_in_T, num_rooks); r_T++){
                    int rooks_left = num_rooks - r_T;
                    long int H_tmp = values[0][rooks_left][fields_number + rooks_left];
                    long int right = single_type_rectangle(i - rooks_left, delta_j, r_T);
                    values[knot_index][num_rooks][fields_number + rooks_left] = H_tmp * right;
                }
            }
        }
        else {
            int maximum_rooks_in_L = min(old_i, delta_j);
            int maximum_rooks_in_T = min(delta_i, delta_j);
            int maximum_rooks_in_W = min(delta_i, old_j);
            for(int num_rooks = min_number_of_rooks(i,j,fields_number); num_rooks <= min(i,j); num_rooks++){
                for(int result = -num_rooks; result <= num_rooks; result++){
                    long int sum = 0;
                    for(int r_L = 0; r_L <= min(maximum_rooks_in_L, num_rooks); r_L++) {
                        for (int r_T = 0; r_T <= min(maximum_rooks_in_T, num_rooks); r_T++) {
                            for (int r_W = 0; r_W <= min(maximum_rooks_in_W, num_rooks); r_W++) {
                                int rooks_left = num_rooks - r_W - r_T - r_L;
                                if(rooks_left >= 0 && min_number_of_rooks(old_i, old_j, fields_number) <= rooks_left && 0 <= fields_number + result - r_W + r_L && result - r_W + r_L < fields_number){
                                    long int H_tmp = values[knot_index - 1][rooks_left][fields_number + result - r_W + r_L];
                                    long int bottom = single_type_rectangle(delta_i, old_j - rooks_left, r_W);
                                    long int corner = single_type_rectangle(delta_i - r_W, delta_j, r_T);
                                    long int right = single_type_rectangle(old_i - rooks_left, delta_j - r_T, r_L);
                                    sum += H_tmp * bottom * corner * right;
                                }
                            }
                        }
                    }
                    values[knot_index][num_rooks][fields_number + result] = sum;
                }
            }
        }
    }
    if(knots.size() > 1){
        i = fields_number;
        j = fields_number;
        vector<int> previous_knot = knots[knots.size() - 2];
        int old_i = previous_knot[0];
        int old_j = previous_knot[1];
        int delta_i = i - old_i;
        int delta_j = j - old_j;
        int maximum_rooks_in_L = min(old_i, delta_j);
        int maximum_rooks_in_T = min(delta_i, delta_j);
        int maximum_rooks_in_W = min(delta_i, old_j);
        int num_rooks = fields_number;
        for(int result = -num_rooks; result <= num_rooks; result++) {
            long int sum = 0;
            for (int r_L = 0; r_L <= min(maximum_rooks_in_L, num_rooks); r_L++) {
                for (int r_T = 0; r_T <= min(maximum_rooks_in_T, num_rooks); r_T++) {
                    for (int r_W = 0; r_W <= min(maximum_rooks_in_W, num_rooks); r_W++) {
                        int rooks_left = num_rooks - r_W - r_T - r_L;
                        if (rooks_left >= 0 && min_number_of_rooks(old_i, old_j, fields_number) <= rooks_left &&
                            0 <= fields_number + result - r_W + r_L && result - r_W + r_L < fields_number) {
                            long int H_tmp = values[knots.size() - 2][rooks_left][fields_number + result - r_W + r_L];
                            long int bottom = single_type_rectangle(delta_i, old_j - rooks_left, r_W);
                            long int corner = single_type_rectangle(delta_i - r_W, delta_j, r_T);
                            long int right = single_type_rectangle(old_i - rooks_left, delta_j - r_T, r_L);
                            sum += H_tmp * bottom * corner * right;
                        }
                    }
                }
            }
            values[knots.size()-1][num_rooks][fields_number + result] = sum;
        }
    }
    return values.back().back();
}

void test_H_chopstick(int A, int n){
    vector<vector<int>> strats = divides_strategies(A,n);
    int errors = 0;
    int tries = 0;
    for (int i = 0; i < strats.size(); i++){
        for (int j = 0; j < strats.size(); j++) {
            tries += 1;
            vector<int> mock_A = strats[i];
            vector<int> mock_B = strats[j];
            vector<long int> h_perms = h_chopstick(mock_A, mock_B);
            vector<long int> h_recc = recurrence_H_chopstick(mock_A, mock_B);
            for(int i = 0; i < h_perms.size(); i++){
                if(h_perms[i] != h_recc[i]) {
                    cout << "error for pair below" << endl;
                    print_vector(mock_A);
                    print_vector(mock_B);
                    errors++;
                    break;
                }
            }
        }
    }
    cout << "Number of errors " << errors <<" on " << tries << " tries." << endl;
}

double symmetrized_payoff_chopstick(vector<int> s_A, vector<int> s_B){
    int fields_number = s_A.size();
    vector<long int> values = recurrence_H_chopstick(s_A, s_B);
    double payoff = 0;
    for (int result = -fields_number; result <= fields_number; result++){
        payoff += values[fields_number + result] * sign(result);
    }
    payoff /= factorial(fields_number);
    return payoff;
}

tuple<vector<vector<double>>, vector<vector<int>>, vector<vector<int>>> payoff_matrix_chopstick_parallel(int A, int B, int n) {
    if(A == B){
        cout << "symmetric" << endl;
        return payoff_matrix_chopstick_symmetric_parallel(A, B, n);
    }
    cout << "unsymmetric" << endl;
    return payoff_matrix_chopstick_unsymmetric_parallel(A, B, n);
}

tuple<vector<vector<double>>, vector<vector<int>>, vector<vector<int>>>  payoff_matrix_chopstick_unsymmetric_parallel(int A, int B, int n){
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
            matrix.at(i).at(j) = symmetrized_payoff_chopstick(A_strategies.at(i), B_strategies.at(j));
        }
    }
    return tie(matrix, A_strategies, B_strategies);
}

tuple<vector<vector<double>>, vector<vector<int>>, vector<vector<int>>>  payoff_matrix_chopstick_symmetric_parallel(int A, int B, int n){
    vector<vector<double>> matrix;
    vector<vector<int>> A_strategies = divides_strategies(A, n);
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
            double payoff =  symmetrized_payoff_chopstick(A_strategies.at(i), A_strategies.at(j));
            matrix.at(i).at(j) = payoff;
            matrix.at(j).at(i) = -payoff;
        }
    }
    return tie(matrix, A_strategies, A_strategies);
}
