#pragma once
class Macierz
{
public:
	double** M;
	int n;

	Macierz(int a1, int a2, int a3, int N);
	Macierz(const Macierz& M);

	~Macierz();
};

