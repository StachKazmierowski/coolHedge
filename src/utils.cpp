#include "utils.h"

int sign(int a){
    if(a > 0)
        return 1;
    if(a == 0)
        return 0;
    return -1;
}

vector<vector<int>> empty_matrix(int size){
    vector<vector<int>> values;
    for(int i = 0 ; i < size; i++){
        vector<int> row;
        for(int j = 0; j < size; j++){
            row.push_back(0);
        }
        values.push_back(row);
    }
    return values;
}

vector<vector<int>> clash_matrix(vector<int> s_A, vector<int> s_B){
    int fields_number = s_A.size();
    vector<vector<int>> clash_matrix;
    clash_matrix = empty_matrix(fields_number);
    for(int i = 0; i < fields_number; i++){
        for(int j = 0; j < fields_number; j++){
            clash_matrix.at(i).at(j) = sign(s_A.at(j) - s_B.at(i));
        }
    }
    return clash_matrix;
}

void print_matrix(vector<vector<int>> matrix){
    cout << "====================" << endl;
    int rows_number = matrix.size();
    int columns_number = matrix.at(0).size();
    for(int i = 0; i < rows_number; i++){
        for(int j = 0; j < columns_number; j++){
            cout << matrix.at(i).at(j) << ",";
        }
        cout << endl;
    }
}

void print_matrix(vector<vector<double>> matrix){
    cout << "====================" << endl;
    int rows_number = matrix.size();
    int columns_number = matrix.at(0).size();
    for(int i = 0; i < rows_number; i++){
        for(int j = 0; j < columns_number; j++){
            cout << matrix.at(i).at(j) << ",";
        }
        cout << endl;
    }
}

void print_vector(vector<int> v){
    cout << "====================" << endl;
    for(int i = 0; i < v.size() ; i++){
        cout << v.at(i) << ",";
    }
    cout << endl;
}

bool check_if_contains(vector<vector<int>> vectors, vector<int> vector){
    if(find(vectors.begin(), vectors.end(), vector) != vectors.end()) {
        return true;
    } else {
        return false;
    }
}

vector<vector<int>> next_divide(vector<vector<int>> divide){
    int fields_number = divide.at(0).size();
    int div_num = divide.size();
    vector<vector<int>> new_divides;
    for(int i =0; i < div_num; i++){
        vector<int> curr_divide = divide.at(i);
        for(int j = 0; j < fields_number; j++){
            if (j == 0 || curr_divide.at(j) < curr_divide.at(j-1)){
                vector<int> new_divide;
                new_divide = curr_divide;
                new_divide[j]++;
                if(!check_if_contains(new_divides, new_divide))
                    new_divides.push_back(new_divide);
            }
        }
    }
    return new_divides;
}

vector<vector<int>> divides_strategies(int resources, int fields_number){
    vector<vector<int>> divides_result;
    vector<int> first_divide;
    for(int i = 0 ; i < fields_number; i++)
        first_divide.push_back(0);
    divides_result.push_back(first_divide);
    for(int i = 0 ; i < resources; i++){
        divides_result = next_divide(divides_result);
    }
    return divides_result;
}

vector<vector<int>> all_permutations(vector<int> s_A){
    reverse(s_A.begin(), s_A.end());
    vector<vector<int>> result;
    do {
        result.push_back(s_A);
    } while (next_permutation(s_A.begin(), s_A.end()));
    return result;
}

int k_W(vector<int> s_A, vector<int> s_B){
    int k_W = 0;
    for (int i = 0; i < s_A.size(); i++){
        if(s_A.at(i) > s_B.at(i))
            k_W++;
    }
    return k_W;
}

int k_L(vector<int> s_A, vector<int> s_B){
    int k_L = 0;
    for (int i = 0; i < s_A.size(); i++){
        if(s_A.at(i) < s_B.at(i))
            k_L++;
    }
    return k_L;
}

int factorial(int n){
    if(n < 0)
        return 0;
    if(n == 0 || n == 1)
        return 1;
    return n* factorial(n-1);
}

int matrix_sum(vector<vector<int>> matrix){
    int sum = 0;
    for(int i = 0; i < matrix.size(); i++){
        for(int j = 0; j < matrix.at(0).size(); j++){
            sum += matrix.at(i).at(j);
        }
    }
    return sum;
}

vector<vector<int>> h(vector<int> s_A, vector<int> s_B){
    int fields_number = s_A.size();
    vector<vector<int>> values;
    values = empty_matrix(fields_number + 1);
    vector<vector<int>> permutations;
    permutations = all_permutations(s_A);
    for (int i = 0; i < permutations.size(); i++){
        vector<int> permutation = permutations.at(i);
        int kW = k_W(permutation, s_B);
        int kL = k_L(permutation, s_B);
        values.at(kW).at(kL)++;
    }
    int scale = int(factorial(fields_number) / matrix_sum(values));
    for(int i = 0; i < values.size(); i++){
        for(int j = 0; j < values.at(0).size(); j++){
            values.at(i).at(j) *= scale;
        }
    }
    return values;
}

vector<int> find_area_vector(vector<vector<int>> clash_mat, int area){
    vector<int> L;
    int fields_number = clash_mat.size();
    for (int i = 0; i < fields_number ; i++){
        int cells_in_L = 0;
        for (int j = 0; j < fields_number ; j++){
            if(clash_mat.at(j).at(i) == area)
                cells_in_L++;
        }
        L.push_back(cells_in_L);
    }
    return L;
}

bool has_single_value(vector<int> v){
    if(set<int>(v.begin(), v.end()).size() == 1)
        return true;
    return false;
}

int number_of_same_values_at_tail(vector<int> v){
    int width = 0;
    for (int i = v.size() - 1; i >= 0; i--){
        if(v.at(i) == v.back())
            width++;
        else
            break;
    }
    return width;
}

int height_to_cutoff(vector<int> L, vector<int> T, int current_height){
    if(L.back() + T.back() < current_height)
        return current_height - (L.back() + T.back());
    if(L.back() == current_height)
        return 0;
    if(T.back() > 0){
        if (has_single_value(T) && has_single_value(L))
            return T.back();
        else if (has_single_value(L) && L.back() == 0 && T.at(0) < current_height)
            return 0;
        else
            return T.back();
    }
    cerr << ("height_to_cutoff FUNCTION FAILED TO RETURN VALUES") << endl;
    return -1;
}

vector<int> add_vectors(vector<int> a, vector<int> b){
    vector<int> result;
    for (int i = 0; i < a.size(); i++){
        result.push_back(a.at(i) + b.at(i));
    }
    return result;
}

int width_to_cutoff(vector<int> L, vector<int> T, int current_height){
    if(L.back() + T.back() < current_height)
        return 0;
    if(L.back() == current_height)
        return number_of_same_values_at_tail(L);
    if(T.back() > 0){
        if (has_single_value(T) && has_single_value(L))
            return 0;
        else if (has_single_value(L) && L.back() == 0 && T.at(0) < current_height)
            return number_of_same_values_at_tail(T);
        else
            return min(number_of_same_values_at_tail(L), number_of_same_values_at_tail(add_vectors(L, T)));
    }
    cerr << ("height_to_cutoff FUNCTION FAILED TO RETURN VALUES") << endl;
    return -1;
}

vector<int> trim_vector(vector<int> v, int number_of_elements){
    vector<int> result;
    for (int i = 0; i < number_of_elements; i++){
        result.push_back(v.at(i));
    }
    return result;
}

vector<vector<int>> find_knots(vector<int> L, vector<int> T){ //TODO następny krok
    int fields_number = L.size();
    int i = fields_number;
    int j = fields_number;
    int current_height = fields_number;
    vector<vector<int>> knots;
    vector<int> knot;
    knot.push_back(i);
    knot.push_back(j);
    knots.push_back(knot);
    while (L.size() > 0){
        int height_to_remove = height_to_cutoff(L, T, current_height);
        int width_to_remove = width_to_cutoff(L, T, current_height);
        i -= height_to_remove;
        j -= width_to_remove;
        current_height -= height_to_remove;
        if(width_to_remove > 0){
            L = trim_vector(L, j);
            T = trim_vector(T, j);
        }
        if(i > 0 && j > 0){
            vector<int> knot;
            knot.push_back(i);
            knot.push_back(j);
            knots.push_back(knot);
        } else {
            break;
        }
    }
    reverse(knots.begin(), knots.end());
    return knots;
}

int newton_symbol(int n, int k){
    if(n < 0 || k < 0 || k > n)
        return 0;
    return int(factorial(n) / (factorial(n - k) * factorial(k)));
}

int single_type_rectangle(int cols_num, int rows_num, int rooks_num){
    if(cols_num < 0 or rows_num < 0 or rooks_num < 0)
        return 0;
    if(rooks_num > cols_num or rooks_num > rows_num)
        return 0;
    return newton_symbol(rows_num, rooks_num) * newton_symbol(cols_num, rooks_num) * factorial(rooks_num);
}

int min_number_of_rooks(int i, int j, int n){
    int delta_i = n - i;
    int delta_j = n - j;
    int r_1 = min(delta_i, j);
    int r_3 = min(i, delta_j);
    int r_2 = min(delta_i - r_1, delta_j - r_3);
    return n - (r_1 + r_2 + r_3);
}

vector<vector<vector<vector<int>>>> empty_H(int knots_number, int fields_number){
    vector<vector<vector<vector<int>>>> values;
    for(int i = 0; i < knots_number; i++){
        vector<vector<vector<int>>> values_0;
        for(int j = 0; j < fields_number + 1; j++){
            vector<vector<int>> values_1;
            for(int j = 0; j < fields_number + 1; j++){
                vector<int> values_2;
                for(int j = 0; j < fields_number + 1; j++){
                    values_2.push_back(0);
                }
                values_1.push_back(values_2);
            }
            values_0.push_back(values_1);
        }
        values.push_back(values_0);
    }
    return values;
}

vector<vector<int>> recurrence_H(vector<int> s_A, vector<int> s_B){
    int fields_number = s_A.size();
    vector<vector<int>> clash_mat = clash_matrix(s_A, s_B);
    vector<int> L = find_area_vector(clash_mat, -1);
    vector<int> T = find_area_vector(clash_mat, 0);
    vector<vector<int>> knots = find_knots(L, T);
    vector<vector<vector<vector<int>>>> values = empty_H(knots.size(), fields_number);
    vector<int> knot = knots.at(0);
    int i = knot.at(0);
    int j = knot.at(1);
    if(L.at(0) > 0){
        for (int num_rooks = 0; num_rooks <= min(i,j); num_rooks++)
            values[0][num_rooks][0][num_rooks] = single_type_rectangle(i, j, num_rooks);
    } else if (T.at(0) > 0){
        for (int num_rooks = 0; num_rooks <= min(i,j); num_rooks++)
            values[0][num_rooks][0][0] = single_type_rectangle(i, j, num_rooks);
    } else {
        for (int num_rooks = 0; num_rooks <= min(i,j); num_rooks++)
            values[0][num_rooks][num_rooks][0] = single_type_rectangle(i, j, num_rooks);
    }
    for(int knot_index = 1; knot_index < knots.size(); knot_index++){
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
                    int H_tmp = values[0][rooks_left][0][rooks_left];
                    int bottom = single_type_rectangle(delta_i, j - rooks_left, r_T);
                    values[knot_index][num_rooks][0][rooks_left] = H_tmp * bottom;
                }
            }
        } else if (knot_index == 1 && T[0] == 0 && T[j-1] > 0 && L[0] == 0 && has_single_value(trim_vector(L, j))){ //W and T stripes
            int maximum_rooks_in_T = min(i, delta_j);
            for(int num_rooks = 0; num_rooks <= min(i,j); num_rooks++){
                for(int r_T = 0; r_T <= min(maximum_rooks_in_T, num_rooks); r_T++){
                    int rooks_left = num_rooks - r_T;
                    int H_tmp = values[0][rooks_left][rooks_left][0];
                    int right = single_type_rectangle(i - rooks_left, delta_j, r_T);
                    values[knot_index][num_rooks][rooks_left][0] = H_tmp * right;
                }
            }
        }
        else {
            int maximum_rooks_in_L = min(old_i, delta_j);
            int maximum_rooks_in_T = min(delta_i, delta_j);
            int maximum_rooks_in_W = min(delta_i, old_j);
            for(int num_rooks = min_number_of_rooks(i,j,fields_number); num_rooks <= min(i,j); num_rooks++){
                for(int k_W = 0; k_W <= num_rooks; k_W++){
                    for(int k_L = 0; k_L <= num_rooks - k_W; k_L++){
                        int sum = 0;
                        for(int r_L = 0; r_L <= min(min(maximum_rooks_in_L, num_rooks), k_L); r_L++) {
                            for (int r_T = 0; r_T <= min(maximum_rooks_in_T, num_rooks); r_T++) {
                                for (int r_W = 0; r_W <= min(min(maximum_rooks_in_W, num_rooks), k_W); r_W++) {
                                    int rooks_left = num_rooks - r_W - r_T - r_L;
                                    if(rooks_left >= 0 && min_number_of_rooks(old_i, old_j, fields_number) <= rooks_left){
                                        int H_tmp = values[knot_index - 1][rooks_left][k_W - r_W][k_L-r_L];
                                        int bottom = single_type_rectangle(delta_i, old_j - rooks_left, r_W);
                                        int corner = single_type_rectangle(delta_i - r_W, delta_j, r_T);
                                        int right = single_type_rectangle(old_i - rooks_left, delta_j - r_T, r_L);
                                        sum += H_tmp * bottom * corner * right;
                                    }
                                }
                            }
                        }
                        values[knot_index][num_rooks][k_W][k_L] = sum;
                    }
                }
            }
        }
    }
    return values.back().back();
}

bool assert_matrices(vector<vector<int>> matrix_1, vector<vector<int>> matrix_2){
    if(matrix_1.size() != matrix_2.size())
        return false;
    for(int i = 0; i < matrix_1.size(); i++){
        if(matrix_1.at(i).size() != matrix_2.at(i).size())
            return false;
    }
    for(int i = 0; i < matrix_1.size(); i++){
        for(int j = 0; j < matrix_1.at(i).size(); j++){
            if(matrix_1.at(i).at(j) != matrix_2.at(i).at(j))
                return false;
        }
    }
    return true;
}

bool assert_matrices(vector<vector<double>> matrix_1, vector<vector<double>> matrix_2){
    if(matrix_1.size() != matrix_2.size())
        return false;
    for(int i = 0; i < matrix_1.size(); i++){
        if(matrix_1.at(i).size() != matrix_2.at(i).size())
            return false;
    }
    for(int i = 0; i < matrix_1.size(); i++){
        for(int j = 0; j < matrix_1.at(i).size(); j++){
            if(matrix_1.at(i).at(j) != matrix_2.at(i).at(j))
                return false;
        }
    }
    return true;
}

void test_H(int A, int n){
    vector<vector<int>> strats = divides_strategies(A,n);
    int errors = 0;
    int tries = 0;
    for (int i = 0; i < strats.size(); i++){
        for (int j = 0; j < strats.size(); j++) {
            tries += 1;
            vector<int> mock_A = strats[i];
            vector<int> mock_B = strats[j];
            vector<vector<int>> h_perms = h(mock_A, mock_B);
            vector<vector<int>> h_recc = recurrence_H(mock_A, mock_B);
            if(!assert_matrices(h_perms, h_recc)) {
                cout << "error for pair below" << endl;
                print_vector(mock_A);
                print_vector(mock_B);
                errors++;
            }
        }
    }
    cout << "Number of errors " << errors <<" on " << tries << " tries." << endl;
}

int blotto(int k_W, int k_L, int n){
    return k_W - k_L;
}

int attack(int k_W, int k_L, int n){
    if(k_W > 0)
        return 1;
    return -1;
}

int majoritarian(int k_W, int k_L, int n){
    return sign(k_W - k_L);
}

int chopstick(int k_W, int k_L, int n){
    if(k_W > n/2)
        return 1;
    if(k_L > n/2)
        return -1;
    return 0;
}


double symmetrized_payoff(vector<int> s_A, vector<int> s_B, int (*aggregation_function)(int, int, int)){
    int fields_number = s_A.size();
    vector<vector<int>> values = recurrence_H(s_A, s_B);
    double result = 0;
    for (int k_W = 0; k_W < values.size(); k_W++){
        for (int k_L = 0; k_L < values.size(); k_L++) {
            result += values[k_W][k_L] * aggregation_function(k_W, k_L, fields_number);
        }
    }
    result /= factorial(fields_number);
    return result;
}

bool is_symmetric(int (*aggregation_function)(int, int, int), int fields_number){
    for (int i = 0; i <= fields_number; i++){
        for (int j = 0; j <= fields_number - i; j++){
            if(aggregation_function(i, j, fields_number) != -aggregation_function(j, i, fields_number)) {
                cout << i << ", " << j << endl;
                cout << aggregation_function(i, j, fields_number) << endl;
                cout << -aggregation_function(j, i, fields_number) << endl;
                return false;
            }
        }
    }
    return true;
}

vector<vector<double>> payoff_matrix(int A, int B, int n, int (*aggregation_function)(int, int, int)){
    vector<vector<double>> matrix;
    vector<vector<int>> A_strategies = divides_strategies(A, n);
    vector<vector<int>> B_strategies = divides_strategies(B, n);
    for(int i = 0 ; i < A_strategies.size(); i++){
        vector<double> row;
        for(int j = 0 ; j < B_strategies.size(); j++){
            row.push_back(symmetrized_payoff(A_strategies.at(i), B_strategies.at(j), aggregation_function));
        }
        matrix.push_back(row);
    }
    return matrix;
}

vector<vector<double>> payoff_matrix_parallel(int A, int B, int n, int (*aggregation_function)(int, int, int)){
    if(is_symmetric(aggregation_function, n)){
        cout << "symmetric" << endl;
        return payoff_matrix_symmetric_parallel(A, B, n, aggregation_function);
    }
    cout << "unsymmetric" << endl;
    return payoff_matrix_unsymmetric_parallel(A, B, n, aggregation_function);
}

vector<vector<double>> payoff_matrix_unsymmetric_parallel(int A, int B, int n, int (*aggregation_function)(int, int, int)){
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
    #pragma omp parallel for default(none) private(i, j) shared(A_strategies, B_strategies, aggregation_function, matrix)
    for(i = 0 ; i < A_strategies.size(); i++){
        for(j = 0 ; j < B_strategies.size(); j++){
            matrix.at(i).at(j) = symmetrized_payoff(A_strategies.at(i), B_strategies.at(j), aggregation_function);
        }
    }
    return matrix;
}

vector<vector<double>> payoff_matrix_symmetric_parallel(int A, int B, int n, int (*aggregation_function)(int, int, int)){
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
#pragma omp parallel for default(none) private(i, j) shared(A_strategies, B_strategies, aggregation_function, matrix)
    for(i = 1 ; i < A_strategies.size(); i++){
        for(j = i + 1 ; j < B_strategies.size(); j++){
            double payoff =  symmetrized_payoff(A_strategies.at(i), B_strategies.at(j), aggregation_function);
            matrix.at(i).at(j) = payoff;
            matrix.at(j).at(i) = -payoff;
        }
    }
    return matrix;
}

