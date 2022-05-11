// Julia Scheffler, 175718

#include <iostream>
#include <math.h>
#include "Macierz.h"
#include "UkladRownan.h"

using namespace std;

int main() {
	double norma;
	cout << "ZADANIE A" << endl;
	// tworzenie uk³adu równañ Ax = b
	int N = 918;
	int e = 7;
	int f = 5;
	// tworzenie macierzy A o wymiarach NxN
	Macierz* Aa = new Macierz(5 + e, -1, -1, N);
	// tworzenie wektora b o d³ugoœci N
	double* ba = new double[N];

	for (int i = 0; i < N; i++) {
		ba[i] = sin(i * (f + 1));
	}

	double* xa = new double[N];

	UkladRownan* ukladRownanA = new UkladRownan(Aa, ba, xa);

	cout << endl << "ZADANIE B" << endl;
	// metoda Jacobiego
	ukladRownanA->metodaJacobiego(&norma);
	cout << "Metoda Jacobiego (N = " << N  << ", a1 = " << 5 + e << ")" << endl;
	cout << "\tLiczba iteracji: " << ukladRownanA->liczbaIteracjiJ << endl;
	cout << "\tCzas: " << ukladRownanA->czasJ << endl;
	cout << "\tNorma z residuum: " << norma << endl << endl;
	// metoda Gaussa-Seidla
	ukladRownanA->metodaGaussaSeidla(&norma);
	cout << "Metoda Gaussa-Seidla (N = " << N << ", a1 = " << 5 + e << ")" << endl;
	cout << "\tLiczba iteracji: " << ukladRownanA->liczbaIteracjiGS << endl;
	cout << "\tCzas: " << ukladRownanA->czasGS << endl;
	cout << "\tNorma z residuum: " << norma << endl << endl;

	cout << endl << "ZADANIE C" << endl;
	// tworzenie macierzy A o wymaiarach NxN
	Macierz* Ac = new Macierz(3, -1, -1, N);
	// tworzenie wektora b o d³ugoœci N
	double* bc = new double[N];

	for (int i = 0; i < N; i++) {
		bc[i] = sin(i * (f + 1));
	}

	double* xc = new double[N];

	UkladRownan* ukladRownanC = new UkladRownan(Ac, bc, xc);
	// sprawdzenie, czy metody iteracyjne zbiegaj¹ siê dla tego przypadku
	ukladRownanC->metodaJacobiego(&norma);
	cout << "Metoda Jacobiego (N = " << N << ", a1 = " << 3 << ")" << endl;
	cout << "\tLiczba iteracji: " << ukladRownanC->liczbaIteracjiJ << endl;
	cout << "\tCzas: " << ukladRownanC->czasJ << endl;
	cout << "\tNorma z residuum: " << norma << endl << endl;
	ukladRownanC->metodaGaussaSeidla(&norma);
	cout << "Metoda Gaussa-Seidla (N = " << N << ", a1 = " << 3 << ")" << endl;
	cout << "\tLiczba iteracji: " << ukladRownanC->liczbaIteracjiGS << endl;
	cout << "\tCzas: " << ukladRownanC->czasGS << endl;
	cout << "\tNorma z residuum: " << norma << endl << endl;

	cout << endl << "ZADANIE D" << endl;
	// metoda faktoryzacji LU
	ukladRownanC->faktoryzacjaLU(&norma);
	// ile wynosi norma z residuum?
	cout << "Metoda faktoryzaji LU (N = " << N << ", a1 = " << 3 << ")" << endl;
	cout << "\tCzas: " << ukladRownanC->czasLU << endl;
	cout << "\tNorma z residuum: " << norma << endl << endl;

	cout << endl << "ZADANIE E" << endl;
	// uk³ad równañ z zadania A
	// N = {100, 500, 1000, 2000, 3000}
	int n[5] = { 100, 500, 1000, 2000, 3000 };
	cout << "(a1 = " << 5 + e << ")" << endl;
	cout << "\tJacobi\tGS\tLU" << endl;
	cout << "-------------------------------" << endl;
	for (int i = 0; i < 5; i++) {
		Macierz* Ae = new Macierz(5 + e, -1, -1, n[i]);
		double* be = new double[n[i]];

		for (int j = 0; j < n[i]; j++) {
			be[j] = sin(j * (f + 1));
		}

		double* xe = new double[n[i]];
		UkladRownan* ukladRownanE = new UkladRownan(Ae, be, xe);
		// metoda Jacobiego
		ukladRownanE->metodaJacobiego(&norma);
		// metoda Gaussa-Seidla
		ukladRownanE->metodaGaussaSeidla(&norma);
		// metoda faktoryzacji LU
		double norma;
		ukladRownanE->faktoryzacjaLU(&norma);
		// porównanie czasu trwania w zale¿noœci od N - wykresy

		cout << n[i] << "||\t" << ukladRownanE->czasJ << "\t" << ukladRownanE->czasGS << "\t" << ukladRownanE->czasLU << endl;
	}

	return 0;
}