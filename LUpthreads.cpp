#include <iostream>
#include <cstdlib>
#include <pthread.h>
#include <chrono>
#include <cmath>

using namespace std;
using namespace chrono;

// Struct to hold parameters for each thread
struct ThreadData {
    double **A, ** L, ** U;
    int n;
    int start_row;
    int end_row;
    int k;
};

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
void *lu_decomposition(void *args) {
    ThreadData *data = (ThreadData *)args;
    // LU decomposition algorithm
    // Implementation of the pseudocode provided in the assignment
    // Parallelized with Pthreads
    // Each thread works on its assigned portion of the matrix
    for(int i = data->start_row; i < data->end_row; ++i){
        for(int j = data->k; j < data->n; ++j){
            data->A[i][j] = data->A[i][j] - data->L[i][data->k]*data->U[data->k][j];
        }
    }
    return NULL;
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
    
    // Allocate memory for the matrix and initialize with random values
    double **A = new double*[n], **L = new double*[n], **U = new double*[n], ** A_old = new double*[n];
    int * pi = new int[n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
        A_old[i] = new double[n];
        for (int j = 0; j < n; j++) {
            // Initialize with random values
            A[i][j] = drand48();
            A_old[i][j] = A[i][j];
        }
    }
    for(int i = 0; i < n; ++i) pi[i] = i;

    for(int i = 0; i < n; ++i){
        L[i] = new double[n], U[i] = new double[n];
        for(int j = 0; j < n; ++j){
            if(i > j) U[i][j] = 0;
            if(i == j) L[i][i] = 1;
            if(j > i) L[i][j] = 0;
        }
    }

    auto start_time = high_resolution_clock::now();
    // Create an array of pthreads and thread data structures
    pthread_t threads[num_threads];
    ThreadData thread_data[num_threads];

    try{
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
                // Divide the work among threads and create them
            int chunk_size = (n - k) / num_threads;
            for (int i = 0; i < num_threads; i++) {
                thread_data[i].A = A;
                thread_data[i].L = L;
                thread_data[i].U = U;
                thread_data[i].n = n;
                thread_data[i].k = k;
                thread_data[i].start_row = k + i * chunk_size;
                thread_data[i].end_row = k + (i + 1) * chunk_size;
                if (i == num_threads - 1) {
                    // Last thread may handle extra rows if n is not divisible by num_threads
                    thread_data[i].end_row = n;
                }
                pthread_create(&threads[i], NULL, lu_decomposition, (void *)&thread_data[i]);
            }

            // Join threads
            for (int i = 0; i < num_threads; i++) {
                pthread_join(threads[i], NULL);
            }
        }
    }catch(exception E){
        cout<<"Error"<<endl;
    }
    auto end_time = high_resolution_clock::now();
    duration<double> timeTaken = end_time - start_time;

        // Calculate and print the execution time
    cout << "Execution time: " << timeTaken.count() << " seconds" << endl;


    // Verification: Calculate the L2,1 norm of the residual matrix
    double residual_norm = calculate_residual_norm(A_old, L, U, n, pi);
    cout << "L2,1 norm of the residual matrix: " << residual_norm << endl;

    
    // Free allocated memory
    freeMemory(A, n), freeMemory(L, n), freeMemory(U, n), freeMemory(A_old, n);
    freeMemory(pi);
    return 0;
}
