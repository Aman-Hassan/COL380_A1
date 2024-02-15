#include <iostream>
#include <cstdlib>
#include <omp.h>
#include <chrono>
#include <cmath>

using namespace std;
using namespace chrono;

void freeMemory(double** matrix, int n)
{                                                   // deleting all variables with significant memory size
    for(int i = 0; i < n; ++i)
        delete[] matrix[i];
    delete[] matrix;
}

void freeMemory(int * vector){
    delete[] vector;
}

// Function to perform LU decomposition with partial pivoting
void lu_decomposition(double **A, double **L, double **U, int *pi, int n, int num_threads) {
    // LU decomposition algorithm
    // Implementation of the pseudocode provided in the assignment
    // Parallelized with OpenMP directives

    for(int k = 0; k < n; ++k){
        double mx = 0;
        int k_dash = k;
        for(int i = k; i < n; ++i){
            if(mx < abs(A[i][k])){
                mx = abs(A[i][k]);
                k_dash = i;
            }
        }
        if(mx == 0) throw -1;
        swap(pi[k], pi[k_dash]);
        swap(A[k], A[k_dash]);
        for(int i = 0; i <= k-1; ++i){
            swap(L[k][i], L[k_dash][i]);
        }
        U[k][k] = A[k][k];
        for(int i = k+1; i < n; ++i){
            L[i][k] = A[i][k]/U[k][k];
            U[k][i] = A[k][i];
        }
        int i;
        // parallelisation part
        #pragma omp parallel for num_threads(num_threads) private(i) shared(A, L, U)
            for(i = k+1; i < n; ++i){
                for(int j = k+1; j < n; ++j){
                    A[i][j] = A[i][j] - L[i][k]*U[k][j];
                }
            }
    }
}

// Function to calculate the L2,1 norm of the residual matrix
double calculate_residual_norm(double **A, double **L, double **U, int n, int *pi) {
    double res = 0.0;
    for(int j = 0; j < n; ++j){
        double col_residue = 0.0;
        for(int i = 0; i < n; ++i){
            double temp = A[pi[i]][j];
            for(int k = 0; k < n; ++k) temp -= (L[i][k])*(U[k][j]);
            col_residue += temp*temp;
        }
        res += sqrt(col_residue);
    }
    return res;
}

int main(int argc, char *argv[]) {
    // Parse command line arguments for matrix size and number of threads
    if (argc != 3) {
        cout << "Usage: " << argv[0] << " <matrix_size> <num_threads>" << endl;
        return 1;
    }
    int n = atoi(argv[1]);
    int num_threads = atoi(argv[2]);
    // Set the number of OpenMP threads
    // omp_set_num_threads(num_threads);
    
    // Allocate memory for the matrix and initialize with random values
    double **A = new double*[n], ** L = new double*[n], ** U = new double*[n], **A_old = new double*[n];
    int *pi = new int[n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
        A_old[i] = new double[n];
        for (int j = 0; j < n; j++) {
            // Initialize with random values
            A[i][j] = drand48();
            A_old[i][j] = A[i][j];
        }
    }
    // cout<<"Matrix A : \n";
    // for(int i  =0 ; i < n; ++i){
    //     for(int j = 0 ; j < n; ++j){
    //         cout<<A[i][j]<<" ";
    //     }
    //     cout<<'\n';
    // }
    // Initialise the permutation vector, L and U matrices as given in pseudocode
    for(int i = 0; i < n; ++i) pi[i] = i;

    for(int i = 0; i < n; ++i){
        L[i] = new double[n], U[i] = new double[n];
        for(int j = 0; j < n; ++j){
            if(i > j) U[i][j] = 0;
            if(i == j) L[i][i] = 1;
            if(j > i) L[i][j] = 0;
        }
    }
    
    // Perform LU decomposition with parallel OpenMP implementation
    // double start_time = omp_get_wtime();
    auto start_time = high_resolution_clock::now();
    try{
        lu_decomposition(A, L, U, pi, n, num_threads);
    }catch(exception E){
        cout<<"Error"<<endl;
    }
    auto end_time = high_resolution_clock::now();
    duration<double> timeTaken = end_time - start_time;
    
    // Verification: Calculate the L2,1 norm of the residual matrix
    // double residual_norm = calculate_residual_norm(A_old, L, U, n, pi);
    // cout << "L2,1 norm of the residual matrix: " << residual_norm << endl;
    
    // Calculate and print the execution time
    cout << "Execution time: " << timeTaken.count() << " seconds" << endl;

    // cout<<"Matrix L : \n";
    // for(int  i =0 ; i  < n; ++i){
    //     for(int j  =0 ; j < n; ++j){
    //         cout<<L[i][j]<<" ";
    //     }
    //     cout<<'\n';
    // }

    // cout<<"Matrix U : \n";
    // for(int i = 0; i < n; ++i){
    //     for(int j = 0; j < n; ++j){
    //         cout<<U[i][j]<<" ";
    //     }
    //     cout<<'\n';
    // }

    // cout<<"Permutation : \n";
    // for(int i = 0; i < n; ++i){
    //     for(int j = 0; j < n; ++j){
    //         if(pi[i] == j) cout<<1<<" ";
    //         else cout<<0<<" ";
    //     }
    //     cout<<'\n';
    // }

    // Free allocated memory
    freeMemory(A, n), freeMemory(L, n), freeMemory(U, n), freeMemory(A_old, n);
    freeMemory(pi);
    return 0;
}
