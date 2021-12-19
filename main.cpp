#include <iostream>
using namespace std;

#include "src/utils.h"

int main() {
    // TODO check_if_symmetric done
    // TODO symmetryczna macierz done
    // TODO równoległe obliczanie (omp parralel?) done
    // TODO zapisywanie macierzy

//    test_H(15,5);
//    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//    vector<vector<double>> mat = payoff_matrix(13,13,10,chopstick);
//    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

    std::chrono::steady_clock::time_point begin_2 = std::chrono::steady_clock::now();
    vector<vector<double>> mat_1 = payoff_matrix_parallel(12,12,10,chopstick);
    std::chrono::steady_clock::time_point end_2 = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end_2 - begin_2).count() << "[ms]" << std::endl;
//    cout << assert_matrices(mat, mat_1) << endl;
//    print_matrix(mat);
//    print_matrix(mat_1);
    return 0;
}
