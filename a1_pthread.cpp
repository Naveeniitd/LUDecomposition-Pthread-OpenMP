
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <pthread.h>

#include <fstream>

using namespace std;
struct ThreadData {
    int startRow, endRow;
    int pivotRow, currentCol;
    double** a;
    double** l;
    double** u;
    int n;
};
void printMatrix(double** a, int n) {
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}

int findPivotRow(double** a, int k, int n) {
    double maxVal = 0.0;
    int pivotRow = k; // Start with the current row as the default pivot
    for (int i = k; i < n; ++i) {
        if (fabs(a[i][k]) > maxVal) {
            maxVal = fabs(a[i][k]);
            pivotRow = i; // Update the pivot row index if a larger value is found
        }
    }
    return pivotRow; // Return the index of the pivot row
}
void updateUMatrixDirectly(double** a, double** l, double** u, int k, int n) {
    for (int j = k; j < n; ++j) {
        u[k][j] = a[k][j]; // Copying the current row of A to the U matrix
    }
    // Update L matrix for the elements below diagonal in column k
    for (int i = k+1; i < n; ++i) {
        l[i][k] = a[i][k] / u[k][k];
    }
}

// Allocate memory for a matrix
double** allocateMatrix(int n) {
    double** matrix = new double*[n];
    for (int i = 0; i < n; ++i) {
        matrix[i] = new double[n];
    }
    return matrix;
    
}

// Free allocated memory for a matrix
void freeMatrix(double** matrix, int n) {
    for (int i = 0; i < n; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

// Initialize matrix with random numbers
void init(double** a, int n) {
    srand(time(nullptr)); // Seed the random number generator
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i][j] = drand48(); // Generate random numbers between 0 and 1
        }
    }
}

void printLMatrix(double** l, int n) {
    cout << "L matrix:" << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            cout << l[i][j] << " ";
        }
        cout << endl;
    }
}

// Print U matrix
void printUMatrix(double** u, int n) {
    cout << "U matrix:" << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j >= i) cout << u[i][j] << " ";
            else cout << "0 ";
        }
        cout << endl;
    }
}
double** MatrixMult(double** A, double** B, int n) {
    double** C = allocateMatrix(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            C[i][j] = 0.0;
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}


double euclideanNorm(double* v, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}


double L21Norm(double** R, int n) {
    double norm = 0.0;
    for (int j = 0; j < n; ++j) {
        double* column = new double[n];
        for (int i = 0; i < n; ++i) {
            column[i] = R[i][j];
        }
        norm += euclideanNorm(column, n);
        delete[] column;
    }
    return norm;
}

double** P_Matrix(const vector<int>& pi, int n) {
    double** P = allocateMatrix(n);
    for (int i = 0; i < n; ++i) {
        // Initialize P[i][j] to 0.0
        fill(P[i], P[i] + n, 0.0);
        P[i][pi[i]] = 1.0;
    }
    return P;
}

double** computePA(double** P, double** A, int n) {
    return MatrixMult(P, A, n); 
}

double** computeLU(double** L, double** U, int n) {
    return MatrixMult(L, U, n); }

double** computeResidual(double** PA, double** LU, int n) {
    double** R = allocateMatrix(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            R[i][j] = PA[i][j] - LU[i][j];
        }
    }
    return R;
}
void* parallelRowSwap(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    double temp;
    for (int i = data->startRow; i < data->endRow; ++i) {
        temp = data->a[data->pivotRow][i];
        data->a[data->pivotRow][i] = data->a[data->currentCol][i];
        data->a[data->currentCol][i] = temp;
    }
    pthread_exit(NULL);
}


void* parallelMatrixUpdate(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    for (int i = data->startRow; i < data->endRow; ++i) {
        if (i > data->currentCol) {
            for (int j = data->currentCol + 1; j < data->n; ++j) {
                data->a[i][j] = data->a[i][j] - data->l[i][data->currentCol] * data->u[data->currentCol][j];
            }
        }
    }
    pthread_exit(NULL);
}

void LU_Decom(double** a, vector<int>& pi, int n, double** l, double** u, int t) {
    pthread_t threads[t];
    ThreadData threadData[t];
    int rowsPerThread, remainder, startRow;

    for (int k = 0; k < n; ++k) {
        // Pivot selection and swapping (sequential)
        int pivotRow = findPivotRow(a, k, n);
        if (pivotRow != k) {
            // Parallel row swap
            rowsPerThread = n / t;
            remainder = n % t;
            startRow = 0;
            for (int i = 0; i < t; ++i) {
                threadData[i] = {startRow, startRow + rowsPerThread + (i < remainder ? 1 : 0), pivotRow, k, a, l, u, n};
                pthread_create(&threads[i], NULL, &parallelRowSwap, (void*)&threadData[i]);
                startRow += rowsPerThread + (i < remainder ? 1 : 0);
            }
            for (int i = 0; i < t; ++i) {
                pthread_join(threads[i], NULL);
            }
        }


        updateUMatrixDirectly(a, l, u, k, n);

        int rowsPerThread = (n - (k + 1)) / t;
        int startRow = k + 1;
        for (int i = 0; i < t; ++i) {
            int endRow = startRow + rowsPerThread + (i < (n - (k + 1)) % t ? 1 : 0);
            threadData[i] = {startRow, endRow, k, k, a, l, u, n};
            pthread_create(&threads[i], NULL, &parallelMatrixUpdate, (void*)&threadData[i]);
            startRow = endRow;
        }
        for (int i = 0; i < t; ++i) {
            pthread_join(threads[i], NULL);
        }

    }
}


int main(int argc, char* argv[]) {
    if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " <matrix size n> <number of threads t>" << std::endl;
            return 1;
        }


    int n = std::atoi(argv[1]);
    int t = std::atoi(argv[2]);

    if (n <= 0 || t <= 0) {
        std::cerr << "Matrix size and number of threads must be positive integers." << std::endl;
        return 1;
    }
    ofstream resultsFile("results.csv", std::ios::app);
    if (!resultsFile.is_open()) {
           cerr << "Failed to open the file for writing.\n";
           return 1;
       }
    
    std::cout << "Matrix size: " << n << ", Number of threads: " << t << std::endl;
    
    double** a = allocateMatrix(n);
    double** l = allocateMatrix(n);
    double** u = allocateMatrix(n);
    
    vector<int> pi(n);
    init(a, n);
    

    
    auto start = chrono::high_resolution_clock::now();
    
    LU_Decom(a, pi, n, l, u, t);


    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    
    cout << "LU decomposition took " << duration.count() << " milliseconds." << endl;
    resultsFile << "PTHREAD" << ", " << n << ", " << (t == 0 ? "NA" : std::to_string(t)) << ", " << duration.count()<< "\n";
    
    resultsFile.close();
    ///printMatrix(a,n);
    ///printLMatrix(l,n);
    ///printUMatrix(u,n);
    
//    // Compute permutation matrix P from permutation vector pi
    double** p = P_Matrix(pi, n);
    ///printMatrix(p,n);
//
//    // Compute matrices PA and LU
    double** PA = computePA(p, a, n);
   double** LU = computeLU(l, u, n);
//
//    // Compute the residual R = PA - LU
    double** R = computeResidual(PA, LU, n);
//
//    // Compute the L2,1 norm of the residual
    double L21norm = L21Norm(R, n); // Assuming L21Norm is adapted for double**
//
    cout << "L2,1 norm of the residual: " << L21norm << endl;
    freeMatrix(a, n);
    freeMatrix(l, n);
    freeMatrix(u, n);
    freeMatrix(p, n);
    freeMatrix(PA, n);
    freeMatrix(LU, n);
    freeMatrix(R, n);
    
    return 0;
}


