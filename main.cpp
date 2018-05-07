#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#define N 995
#define A1 9
#define A1C 3
#define A23 -1
#define STOP 100

using namespace std;

void jacobi(double** A, double* b, double* x0, double* x1, double* residuum, int& iterationCounter) {
	for (int i = 0; i < N; i++) {
		double sum1 = 0, sum2 = 0;
		for (int j = 0; j <= i - 1; j++) sum1 += A[i][j] * x0[j];
		for (int j = i + 1; j < N; j++) sum2 += A[i][j] * x0[j];
		x1[i] = (b[i] - sum1 - sum2) / A[i][i];
	}
	for (int i = 0; i < N; i++) {
		double temp = 0;
		for (int j = 0; j < N; j++) temp += A[i][j] * x1[j];
		residuum[i] = temp - b[i];
	}
	iterationCounter++;
}

void gauss(double** A, double* b, double* x0, double* x1, double* residuum, int& iterationCounter) {
	for (int i = 0; i < N; i++) {
		double sum1 = 0, sum2 = 0;
		for (int j = 0; j <= i - 1; j++) sum1 += A[i][j] * x1[j];
		for (int j = i + 1; j < N; j++) sum2 += A[i][j] * x0[j];
		x1[i] = (b[i] - sum1 - sum2) / A[i][i];
	}
	for (int i = 0; i < N; i++) {
		double temp = 0;
		for (int j = 0; j < N; j++) temp += A[i][j] * x1[j];
		residuum[i] = temp - b[i];
	}
	iterationCounter++;
}

void luFactor(double** U, double** L) {
	for(int k = 0; k < N - 1; k++)
		for (int j = k + 1; j < N; j++) {
			L[j][k] = U[j][k] / U[k][k];
			for (int m = k; m < N; m++) U[j][m] = U[j][m] - L[j][k] * U[k][m];
		}
}

void stepBack(double** L, double* b, double* y) {
	for (int i = 0; i < N; i++) {
		y[i] = b[i];
		for (int j = 0; j < i; j++) y[i] -= L[i][j] * y[j];
		y[i] /= L[i][i];
	}
}

void stepForward(double** U, double* x, double* y) {
	for (int i = N - 1; i >= 0; i--) {
		x[i] = y[i];
		for (int j = N - 1; j > i; j--) x[i] -= U[i][j] * x[j];
		x[i] /= U[i][i];
	}
}

void reset(double* x0, double* x1, double* residuum, double& norm) {
	for (int i = 0; i < N; i++) {
		x0[i] = 0;
		x1[i] = 0;
		residuum[i] = 0;
		norm = 1;
	}
}

int main() {
	//Zad A, B
	double** A = new double*[N];
	for (int i = 0; i < N; i++) A[i] = new double[N];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++) A[i][j] = 0.0;

	static double b[N], x0[N], x1[N], residuum[N], norm = 1; //x0 macierz xow poprzednich, x1 macierz xow nowych
	for (int i = 0; i < N; i++) {
		b[i] = sin(i * 10);
		for (int j = 0; j < N; j++) {
			if (i == j) A[i][j] = A1;
			else if (i == j + 1 || i == j - 1) A[i][j] = A23;
			else if (i == j + 2 || i == j - 2) A[i][j] = A23;
		}
	}
	int jacobiIters = 0, gaussIters = 0;
	auto start = chrono::system_clock::now();
	while (norm > pow(10, -9) && norm != 0.0) {
		double result = 0;
		jacobi(A, b, x0, x1, residuum, jacobiIters);
		for (int i = 0; i < N; i++) x0[i] = x1[i];
		for (int i = 0; i < N; i++) result += residuum[i] * residuum[i];
		norm = sqrt(result);
	}
	auto end = chrono::system_clock::now();
	cout << "Jacobi iterates: " << jacobiIters << " norm(residuum): " << norm;
	cout << "Time jacobi: " << ((chrono::duration<double>)(end - start)).count() << " s\n";
	reset(x0, x1, residuum, norm);
	start = chrono::system_clock::now();
	while (norm > pow(10, -9) && norm != 0.0) {
		double result = 0;
		gauss(A, b, x0, x1, residuum, gaussIters);
		for (int i = 0; i < N; i++) x0[i] = x1[i];
		for (int i = 0; i < N; i++) result += residuum[i] * residuum[i];
		norm = sqrt(result);
	}
	end = chrono::system_clock::now();
	cout << "\nGauss iterates: " << gaussIters << " norm(residuum): " << norm;;
	cout << "Time gauss: " << ((chrono::duration<double>)(end - start)).count() << " s\n";


	//zad C
	ofstream file;
	file.open("norma.txt");
	file << "Jacobi:\n";
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) A[i][j] = A1C;
			else if (i == j + 1 || i == j - 1) A[i][j] = A23;
			else if (i == j + 2 || i == j - 2) A[i][j] = A23;
		}
	}
	reset(x0, x1, residuum, norm);
	int jacobiIters2 = 0, gaussIters2 = 0, iter = 0;
	double normJ[STOP], normG[STOP];
	while (norm > pow(10, -9) && norm != 0.0 && iter != STOP) {
		double result = 0;
		jacobi(A, b, x0, x1, residuum, jacobiIters2);
		for (int i = 0; i < N; i++) x0[i] = x1[i];
		for (int i = 0; i < N; i++) result += residuum[i] * residuum[i];
		normJ[iter] = sqrt(result);
		file << "(" << iter+1 << ", " << normJ[iter] << ")\n";
		iter++;
	}
	reset(x0, x1, residuum, norm);
	iter = 0;
	file << "Gauss: \n";
	while (norm > pow(10, -9) && norm != 0.0 && iter != STOP) {
		double result = 0;
		gauss(A, b, x0, x1, residuum, gaussIters2);
		for (int i = 0; i < N; i++) x0[i] = x1[i];
		for (int i = 0; i < N; i++) result += residuum[i] * residuum[i];
		normG[iter] = sqrt(result);
		file << "(" << iter+1 << ", " << normG[iter] << ")\n";
		iter++;
	}


	//zad D
	reset(x0, x1, residuum, norm);
	static double y[N], x[N];
	double** L = new double*[N];
	double** U = new double*[N];
	for (int i = 0; i < N; i++) {
		L[i] = new double[N];
		U[i] = new double[N];
	}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++) {
			if (i == j) L[i][j] = 1.0;
			else L[i][j] = 0.0;
			U[i][j] = A[i][j];
		}
	start = chrono::system_clock::now();
	luFactor(U, L);
	double result = 0;
	stepBack(L, b, y);
	stepForward(U, x, y);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			residuum[i] += A[i][j] * x[j];
		}
		residuum[i] -= b[i];
		result += pow(residuum[i], 2);
	}
	norm = sqrt(result);
	end = chrono::system_clock::now();
	cout << "\nnorm: " << norm << '\n';
	cout << "Time LU: " << ((chrono::duration<double>)(end - start)).count() << " s\n";
	file.close();
	return 0;
}