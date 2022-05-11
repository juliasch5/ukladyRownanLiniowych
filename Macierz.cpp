#include "Macierz.h"

Macierz::Macierz(int a1, int a2, int a3, int N) {
    this->n = N;
    this->M = new double*[n];
    for (int i = 0; i < n; i++) {
        this->M[i] = new double[n];
    }
    for (int i = 0; i < this->n; i++) {
        for (int j = 0; j < this->n; j++) {
            if (i == j) {
                M[i][j] = a1;
            }
            else if (i == j + 1 || j == i + 1) {
                M[i][j] = a2;
            }
            else if (i == j + 2 || j == i + 2) {
                M[i][j] = a3;
            }
            else {
                M[i][j] = 0;
            }
        }
    }
}

Macierz::Macierz(const Macierz& M) {
    this->n = M.n;
    this->M = new double*[this->n];
    for (int i = 0; i < n; i++) {
        this->M[i] = new double[n];
    }
    for (int i = 0; i < this->n; i++) {
        for (int j = 0; j < this->n; j++) {
            this->M[i][j] = M.M[i][j];
        }
    }
}

Macierz::~Macierz() {
    for (int i = 0; i < this->n; i++) {
        delete[] this->M[i];
    }
    delete[] this->M;
}
