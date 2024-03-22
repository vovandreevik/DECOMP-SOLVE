#include <iostream>
#include "Forsythe.h"
#include <vector>
#include <iomanip>

void matrixGeneration(double** matrix, int n) {
    for (int i = 0; i < n; i++) {
        matrix[i] = new double[n];
        for (int j = 0; j < n; j++) {
            matrix[i][j] = 1.0 / ((1 + i) + (1 + j) - 1);
        }
    }
}

void matrixTransposition(double** matrix, double** transposedMatrix, int n) {
    for (int i = 0; i < n; i++) {
        transposedMatrix[i] = new double[n];
        for (int j = 0; j < n; j++) {
            transposedMatrix[i][j] = 0.0;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            transposedMatrix[j][i] = matrix[i][j];
        }
    }
}

void matrixMultiplication(double** A, double** B, double** result, int n) {
    for (int i = 0; i < n; i++) {
        result[i] = new double[n];
        for (int j = 0; j < n; j++) {
            result[i][j] = 0;
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void matrixPrint(double** matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << std::setw(10) << std::setprecision(7) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void matrixClear(double** matrix, int n) {
    for (int i = 0; i < n; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

void vectorGeneration(double* vector, int n) {
    for (int i = 0; i < n; i++) {
        vector[i] = 0.0;
    }

    for (int i = 0; i < n; i++) {
        for (int k = 1; k <= n; k++)
        vector[i] += 1.0 / ((i + 1) + k - 1);
    }
}

void matrixVectorMultiplication(double** matrix, double* vector, double* result, int n) {
    for (int i = 0; i < n; i++) {
        result[i] = 0.0;
        for (int j = 0; j < n; j++) {
            result[i] += matrix[i][j] * vector[j];
        }
    }
}

void vectorPrint(double* vector, int n) {
    for (int i = 0; i < n; i++) {
        std::cout << std::setw(10) << std::setprecision(7) << vector[i] << " ";
    }
    std::cout << std::endl;
}

void vectorClear(double* vector) {
    delete[] vector;
}


double sigmaCalculation(double* vectorA, double* vectorB, int n) {
    double temp = 0.0, sigma = 0.0;
    for (int i = 0; i < n; i++) {
        temp += vectorA[i] * vectorA[i];
        sigma += pow(vectorA[i] - vectorB[i], 2);
    }
    return sqrt(sigma / temp);
}


int main() {
    int start = 4, end = 12, step = 2;
    const int COUNT_OF_STEPS = 5;

    double steps[COUNT_OF_STEPS] = { 0 };
    double condsA[COUNT_OF_STEPS] = { 0 };
    double condsB[COUNT_OF_STEPS] = { 0 }; 
    double sigmas[COUNT_OF_STEPS] = { 0 };
    int counter = 0;

    for (int n = start; n <= end; n += step) {
        std::cout << "N = " << n << ":\n";
        steps[counter] = n;

        std::cout << "\nMatrix 1:\n";
        double** matrixA = new double* [n];
        matrixGeneration(matrixA, n);
        matrixPrint(matrixA, n);

        std::cout << "\nVector 1:\n";
        double* vectorA = new double[n];
        vectorGeneration(vectorA, n);
        vectorPrint(vectorA, n);


        //creation matrixB, vectorB
        double** transposedMatrix = new double* [n];
        matrixTransposition(matrixA, transposedMatrix, n);

        double** matrixB = new double* [n];
        matrixMultiplication(transposedMatrix, matrixA, matrixB, n);

        double* vectorB = new double[n];
        matrixVectorMultiplication(transposedMatrix, vectorA, vectorB, n);


        //Decomp & Solve matrixA
        double condA = 0.0;
        std::vector<int> ipvtA(n);
        int k = 0;
        std::vector<double> newA(n * n);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++, k++)
                newA[k] = matrixA[i][j];
        Decomp(n, newA.data(), &condA, ipvtA.data());
        Solve(n, newA.data(), vectorA, ipvtA.data());

        std::cout << "\nResult 1:\n";
        vectorPrint(vectorA, n);

        std::cout << "\nCond 1 = " << std::setprecision(10) << condA << "\n";
        condsA[counter] = condA;
   

        std::cout << "\n\nMatrix 2\n";
        matrixPrint(matrixB, n);

        std::cout << "\nVector 2\n";
        vectorPrint(vectorB, n);

        //Decomp & Solve matrixB
        double condB = 0.0;
        std::vector<int> ipvtB(n);
        k = 0;
        std::vector<double> newB(n * n);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++, k++)
                newB[k] = matrixB[i][j];

        Decomp(n, newB.data(), &condB, ipvtB.data());
        Solve(n, newB.data(), vectorB, ipvtB.data());

        std::cout << "\nResult 2:\n";
        vectorPrint(vectorB, n);

        std::cout << "\nCond 2 = " << std::setprecision(10) << condB << "\n";
        condsB[counter] = condB;

        //Sigma
        double sigma = sigmaCalculation(vectorA, vectorB, n);
        std::cout << "\n\nSigma = " << sigma << '\n';
        sigmas[counter] = sigma;

        counter++;

        std::cout << "\n--------------------------------------------------------------------------------------------------------\n\n";
        matrixClear(matrixA, n);
        matrixClear(transposedMatrix, n);
        matrixClear(matrixB, n);
        vectorClear(vectorA);
        vectorClear(vectorB);
    }

    std::cout << "\nRESULTS: \n";
    std::cout << std::setw(2) << "N"
        << std::setw(11) << "cond 1"
        << std::setw(20) << "cond 2"
        << std::setw(19) << "sigma"
              << std::endl;

    for (int i = 0; i < COUNT_OF_STEPS; i++) {
        std::cout << std::setw(2) << steps[i]
            << std::setw(20) << std::setprecision(10) << condsA[i]
            << std::setw(20) << std::setprecision(10) << condsB[i]
            << std::setw(20) << std::setprecision(10) << sigmas[i]
            << std::endl;
    }

    return 0;
}