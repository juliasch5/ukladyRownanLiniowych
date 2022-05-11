#include <iostream>
#include <math.h>
#include <ctime>
#include <stdio.h>
#include "UkladRownan.h"

UkladRownan::UkladRownan(Macierz* A, double* b, double* x) {
	this->A = A;
	this->b = b;
	this->x = x;
	this->liczbaIteracjiJ = 0;
	this->liczbaIteracjiGS = 0;
	this->czasJ = 0.0;
	this->czasGS = 0.0;
	this->czasLU = 0.0;
	this->normaResiduum = 1e-9;
}

void UkladRownan::metodaJacobiego(double* norma) {
	int N = A->n;
	double* xPoprzedni = new double[N];
	for (int i = 0; i < N; ++i) {
		xPoprzedni[i] = 1;
	}
	double* res = new double[N];
	double normaRes;

	std::clock_t start;
	start = std::clock();

	while (true) {
		liczbaIteracjiJ += 1;
		for (int i = 0; i < N; ++i) {
			double suma1 = 0;
			for (int j = 0; j < i; ++j) {
				suma1 += A->M[i][j] * xPoprzedni[j];
			}
			double suma2 = 0;
			for (int j = i + 1; j < N; ++j) {
				suma2 += A->M[i][j] * xPoprzedni[j];
			}
			x[i] = (b[i] - suma1 - suma2) / A->M[i][i];
		}
		for (int i = 0; i < N; ++i) {
			xPoprzedni[i] = x[i];
		}

		residuum(res);
		normaRes = normaZResiduum(res);
		if (normaRes < normaResiduum || std::isinf(normaRes)) {
			break;
		}
	}
	czasJ = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	*norma = normaRes;
}

void UkladRownan::metodaGaussaSeidla(double* norma) {
	int N = A->n;
	double* xPoprzedni = new double[N];
	for (int i = 0; i < N; ++i) {
		xPoprzedni[i] = 1;
	}
	double* res = new double[N];
	double normaRes;

	std::clock_t start;
	start = std::clock();

	while (true) {
		liczbaIteracjiGS += 1;
		for (int i = 0; i < N; ++i) {
			double suma1 = 0;
			for (int j = 0; j < i; ++j) {
				suma1 += A->M[i][j] * x[j];
			}
			double suma2 = 0;
			for (int j = i + 1; j < N; ++j) {
				suma2 += A->M[i][j] * xPoprzedni[j];
			}
			x[i] = (b[i] - suma1 - suma2) / A->M[i][i];
		}
		for (int i = 0; i < N; ++i) {
			xPoprzedni[i] = x[i];
		}

		residuum(res);
		normaRes = normaZResiduum(res);
		if (normaRes < normaResiduum || std::isinf(normaRes)) {
			break;
		}
	}
	czasGS = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	*norma = normaRes;
}

void UkladRownan::faktoryzacjaLU(double* norma) {
	int N = A->n;

	std::clock_t start;
	start = std::clock();

	Macierz* L = new Macierz(1, 0, 0, N);
	Macierz* U = new Macierz(*A);
	for (int i = 0; i < N - 1; ++i) {
		for (int j = i + 1; j < N; ++j) {
			L->M[j][i] = U->M[j][i] / U->M[i][i];
			for (int k = i; k < N; ++k) {
				U->M[j][k] -= L->M[j][i] * U->M[i][k];
			}
		}
	}
	double* y = new double[N];
	// L * y = b
	for (int i = 0; i < N; ++i) {
		double suma = 0;
		for (int j = 0; j < i; ++j) {
			suma += L->M[i][j] * y[j];
		}
		y[i] = (b[i] - suma) / L->M[i][i];
	}
	// U * x = y
	for (int i = N - 1; i >= 0; --i) {
		double suma = 0;
		for (int j = i + 1; j < N; ++j) {
			suma += U->M[i][j] * x[j];
		}
		x[i] = (y[i] - suma) / U->M[i][i];
	}
	double* res = new double[N];
	residuum(res);
	*norma = normaZResiduum(res);

	czasLU = (std::clock() - start) / (double)CLOCKS_PER_SEC;
}

// res = Ax - b
void UkladRownan::residuum(double*& r) {
	int N = this->A->n;
	for (int i = 0; i < N; ++i) {
		double suma = 0;
		for (int j = 0; j < N; ++j) {
			suma += this->A->M[i][j] * x[j];
		}
		r[i] = suma;
	}
	for (int i = 0; i < N; ++i) {
		r[i] -= b[i];
	}
}

double UkladRownan::normaZResiduum(double* r) {
	int N = A->n;
	double suma = 0;
	for (int i = 0; i < N; ++i) {
		suma += pow(r[i], 2);
	}
	return sqrt(suma);
}

UkladRownan::~UkladRownan() {
	delete[] A;
	delete[] b;
	delete[] x;
}