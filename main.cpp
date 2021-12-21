#include <iostream>
using namespace std;

#include "src/utils.h"
#include "src/chopstick.h"
#include "src/files_utils.h"

int main(int argc, char *argv[]) {
    // TODO różne rodzaje macieży obliczać na raz
    // TODO jakie omp parralel jest najlepsze
    int A = atoi(argv[1]);
    int B = atoi(argv[2]);
    int n = atoi(argv[3]);
    find_and_save_chopstick_matrix(A, B, n);
    return 0;
}
