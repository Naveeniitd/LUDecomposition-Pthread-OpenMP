#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <omp.h>
#include <fstream>

using namespace std;



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
void printMatrix(double** a, int n) {
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << a[i][j] << " ";
        }
        cout << endl;
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
    return MatrixMult(P, A, n); // Assuming MatrixMult is adapted for double**
}

double** computeLU(double** L, double** U, int n) {
    return MatrixMult(L, U, n); // Assuming MatrixMult is adapted for double**
}

double** computeResidual(double** PA, double** LU, int n) {
    double** R = allocateMatrix(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            R[i][j] = PA[i][j] - LU[i][j];
        }
    }
    return R;
}

void LU_Decom(double** a, vector<int>& pi, int n, double** l, double** u) {
    // Initialize π, l, and u
    for (int i = 0; i < n; ++i) {
        pi[i] = i;
        l[i][i] = 1.0;
        for (int j = 0; j < i; ++j) {
            u[i][j] = 0.0;
            l[j][i] = 0.0;
        }
    }

    for (int k = 0; k < n; ++k) {
        double max = 0.0;
        int k_prime = -1;

        // Find the pivot element (Sequential)
        for (int i = k; i < n; ++i) {
            if (max < std::abs(a[i][k])) {
                max = std::abs(a[i][k]);
                k_prime = i;
            }
        }

        // Check for a singular matrix
        if (max == 0.0) {
            std::cerr << "Error: Singular matrix" << std::endl;
            exit(EXIT_FAILURE);
        }

        // Perform pivoting (Sequential parts)
        std::swap(pi[k], pi[k_prime]);
        std::swap(a[k], a[k_prime]);
        for (int i = 0; i < k; ++i) {
            std::swap(l[k][i], l[k_prime][i]);
        }

        u[k][k] = a[k][k];

        // Compute l and u elements (Parallelizable)
        #pragma omp parallel for
        for (int i = k + 1; i < n; ++i) {
            l[i][k] = a[i][k] / u[k][k];
            u[k][i] = a[k][i];
            for (int j = k + 1; j < n; ++j) {
                #pragma omp atomic
                a[i][j] -= l[i][k] * u[k][j];
            }
        }
    }
}
int main(int argc, char* argv[]) {
    if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " <matrix size n> <number of threads t>" << std::endl;
            return 1;
        }

    // Convert arguments from strings to integers
    int n = std::atoi(argv[1]);
    int t = std::atoi(argv[2]);

    if (n <= 0 || t <= 0) {
        std::cerr << "Matrix size and number of threads must be positive integers." << std::endl;
        return 1;
    }
    omp_set_num_threads(t);
    ofstream resultsFile("results.csv", std::ios::app);
    if (!resultsFile.is_open()) {
           cerr << "Failed to open the file for writing.\n";
           return 1;
       }
    

   
    // Now you can use n and t in your program
    std::cout << "Matrix size: " << n << ", Number of threads: " << t << std::endl;
    
    double** a = allocateMatrix(n);
    double** l = allocateMatrix(n);
    double** u = allocateMatrix(n);
    
    vector<int> pi(n);
    init(a, n);
    
    
    
    auto start = chrono::high_resolution_clock::now();
    LU_Decom(a, pi, n, l, u);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    
    cout << "LU decomposition took " << duration.count() << " milliseconds." << endl;
    
    resultsFile << "OPENMP" << ", " << n << ", " << (t == 0 ? "NA" : std::to_string(t)) << ", " << duration.count()<< "\n";
    resultsFile.close();
    ///printMatrix(a,n);
    ///printLMatrix(l,n);
   ///printUMatrix(u,n);
    
    // Compute permutation matrix P from permutation vector pi
    double** p = P_Matrix(pi, n);
    ///printMatrix(p,n);
//    // Compute matrices PA and LU
    double** PA = computePA(p, a, n);
    double** LU = computeLU(l, u, n);
//
//    // Compute the residual R = PA - LU
    double** R = computeResidual(PA, LU, n);
//
//    // Compute the L2,1 norm of the residual
    double L21norm = L21Norm(R, n); // Assuming L21Norm is adapted for double**

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


