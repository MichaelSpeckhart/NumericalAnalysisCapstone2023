#include <vector>
#include <cstdlib>

#include "tbb/tbb.h"

using namespace std;

/**
 * @brief LU Factorization of a square matrix
 * 
 * @param A 
 * @return pair<vector<vector<double>>,vector<vector<double>>> 
 */
pair<vector<vector<double>>,vector<vector<double>>> lu_factorization_parallel(const vector<vector<double>>& A) {
    if (A.size() != A[0].size()) {
        throw invalid_argument("Error: Matrix must be square nxn");
    }
    const int n = static_cast<int>(A.size());
    vector<vector<double>> L(n, vector<double>(n, 0.0)); 
    vector<vector<double>> U = A; 

    tbb::parallel_for(tbb::blocked_range<int>(0, n), 
        [&] (tbb::blocked_range<int>& range) {
            for (int i = range.begin(); i < range.end(); ++i) {
              //Assign the diagonal of the lower triangular matrix to 1
              L[i][i] = 1.0; 
              for (int j = i+1; j < n; ++j) {
                //Find the factor and stored in the lower triangular matrix
                double factor = U[j][i] / U[i][i]; 
                L[j][i] = factor; 
                for (int k = i; k < n; ++k) {
                  //
                  U[j][k] -= factor * U[i][k]; 
                }
              }
            }
    });
    return std::make_pair(L, U);
}

/**
 * @brief 
 * 
 * @param A 
 * @return true 
 * @return false 
 */
bool gaussian_elimination(std::vector<std::vector<double> > &A) {
    if (A.size() != A[0].size())
    {
        throw std::invalid_argument("The matrix must be square.");
    }
     // Iterate over each row in the matrix
    double pivot;
    for(size_t i = 0; i < A.size() - 1; i++){
        // Pivot will be the diagonal
        pivot = A[i][i];
        //stopwatch.elapsed();
        // Iterate of the remaining row elements
        tbb::parallel_for( tbb::blocked_range<size_t>(i+1, A[0].size()), [&](tbb::blocked_range<size_t> r){
            for(size_t j = r.begin(); j < r.end(); j++){
                A[i][j] /= pivot;
            }
        });

        // Do direct assignment for trivial case (self-divide)
        A[i][i] = 1.0;

        // Eliminate ith element from the jth row
            tbb::parallel_for( tbb::blocked_range<size_t>(i+1, A.size()), [&](tbb::blocked_range<size_t> r){
                float scale;
                for(size_t j = r.begin(); j < r.end(); j++){
                    // Factor we will use to scale subtraction by
                    scale = A[j][i];

                    // Iterate over the remaining columns
                    for(size_t k = i + 1; k < A.size(); k++){
                        A[j][k] -= A[i][k] * scale;
                    }

                    // Do direct assignment for trivial case (eliminate position)
                    A[j][i] = 0;
                }
            });
    }
    A[A.size()-1][A[0].size()-1] = 1;

    return true;
}

/**
 * @brief As a prerequisite for the Jacobi Method, the matrix must be diagonally dominant,
 * meaning that the elements on the diagonal indices of the matrix are greater or equal to
 * the sum of the rest of the elements in that row.
 * 
 * @param denseMatrix 
 * @return true If is diagonally dominant
 * @return false Otherwise
 */
bool diagonally_dominant(std::vector<std::vector<double>> denseMatrix) {
    for (size_t i = 0; i < denseMatrix.size(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < denseMatrix[i].size(); ++j) {
            if (j != i) {
                sum += std::abs(denseMatrix[i][j]);
            }
        }
        if (denseMatrix[i][i] < sum) {
            return false;
        }
    }
    
    return true;
}

/**
 * @brief Jacobi Method parallel that uses a parallel for to map each thread to a row and reduces the summation and product
 * 
 * @param A 
 * @param B 
 * @param maxIterations 
 * @return std::vector<double> 
 */
std::vector<double> jacobi_method_parallel(vector<vector<double>> A, vector<double> B, const int maxIterations) {
    if (diagonally_dominant(A) == false) {
        throw std::invalid_argument("Input matrix is not diagonally dominant");
    }
    //Store the coefficients in two vectors, the updated values and the approximated values
    std::vector<double> xValues(B.size(), 0.0);
    std::vector<double> approxValues(B.size(), 0.0);
    atomic<int> iterations = 0;
    const size_t n = static_cast<size_t>(A.size());
    while (iterations < maxIterations) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
            [&] (tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i < range.end(); ++i) {
                    const size_t m = static_cast<size_t>(A[i].size());
                    double sum = 0.0;
                    for (size_t j = 0; j < m; ++j) {
                        if (j != i) {
                            sum += A[i][j] * xValues[j];
                        }
                    }
                    approxValues[i] = (B[i] - sum) / A[i][i];
                    xValues = approxValues;
                    iterations++;
                }
            });
    
    }
    return approxValues;
}

vector<vector<double> > sum_matrix_parallel(vector<vector<double> > m1,
                                   const vector<vector<double> > m2) {
  //   NB: in the future it may be important to use the compiler defined sizes
  //   for platform portability
  size_t d1 = static_cast<size_t>(m1.size()),
            d2 = static_cast<size_t>(m1[0].size());

  // If the matrices do not have the same dimensions we return a dimension error
  if (d1 != static_cast<size_t>(m2.size()) || d2 != static_cast<size_t>(m2[0].size()))
    throw invalid_argument("Matrices do not have same dimensions");

  tbb::parallel_for(tbb::blocked_range<size_t>(0, d1),
    [&] (tbb::blocked_range<size_t>& range) {
      for (size_t i = range.begin(); i < range.end(); ++i) {
        for (size_t j = 0; j < d2; ++j) {
          m1[i][j] -= m2[i][j];
        }
      }
    });

  return m1;
}