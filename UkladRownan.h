#pragma once
#include "Macierz.h"

class UkladRownan
{
public:
	Macierz* A;
	double* b;
	double* x;
	int liczbaIteracjiJ;
	int liczbaIteracjiGS;
	double czasJ;
	double czasGS;
	double czasLU;
	double normaResiduum;

	UkladRownan(Macierz* A, double* b, double* x);

	void metodaJacobiego(double* norma);
	void metodaGaussaSeidla(double* norma);
	void faktoryzacjaLU(double *norma);

	void residuum(double*& r);
	double normaZResiduum(double* r);

	~UkladRownan();
};

