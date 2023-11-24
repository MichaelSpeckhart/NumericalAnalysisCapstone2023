// functions.cc
// Contains all of the functions for performing operations on matrices.

// Potential improvement, implement the functions using iterators instead of
// counter based loops.
// See: algorithm.h by GCC

#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <climits>
#include <string>
#include <vector>
// #include <omp.h>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <tuple>

typedef std::vector<std::vector<double>> matrix;

using namespace std;

/// error with dimensions in a matrix
vector<vector<double>> d_err{{INT_MAX, INT_MAX, INT_MIN},
                             {INT_MIN, INT_MIN, INT_MAX}};
/// arithmetic error
vector<vector<double>> a_err{{INT_MIN, INT_MIN, INT_MAX},
                             {INT_MAX, INT_MAX, INT_MIN}};

/// @brief Reads a file and transforms the bytes into a 2d vector (a.k.a
/// matrix).
/// @param filename The name of the file with the data.
/// @return 2d vector containing the data.
vector<vector<vector<double>>> read_file(char *filename)
{
    try
    {
        fstream file(filename);
        if (file.fail())
        {
            cout << "error (" << errno << "): failed to open file '" << filename
                 << "'" << endl;
            return vector<vector<vector<double>>>();
        }

        string row;
        vector<string> data;
        while (getline(file, row))
            data.push_back(row);

        string s_temp;
        vector<string> f_data;

        for (auto i : data)
        {
            stringstream stream(i);
            while (getline(stream, s_temp, ','))
                f_data.push_back(s_temp);
        }

        int num;

        int d1, d2;

        vector<vector<double>> temp;
        vector<vector<vector<double>>> p_matrices;

        data = f_data;

        vector<double> line;

        num = (int)stof(data[0]);
        data.erase(data.begin(), data.begin() + 1);

        for (int i = 0; i < num; ++i)
        {
            d1 = (int)stof(data[0]);
            d2 = (int)stof(data[1]);
            data.erase(data.begin(), data.begin() + 2);

            for (int j = 0; j < d1; ++j)
            {
                for (int k = 0; k < d2; ++k)
                {
                    double current = (double)stof(data[k]);
                    line.push_back(current);
                }
                temp.push_back(line);
                line.clear();
                data.erase(data.begin(), data.begin() + d2);
            }
            p_matrices.push_back(temp);
            temp.clear();
        }
        return p_matrices;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return {{}};
    }
}

/// @brief Generates an identity matrix.
/// @param size The size of the matrix.
/// @return The identity matrix.
vector<vector<double>> identity_matrix(size_t size)
{
    vector<vector<double>> matrix(size, vector<double>(size, 0));
    for (size_t i = 0; i < size; ++i)
        matrix[i][i] = 1;
    return matrix;
}

/// @brief Adds to matrices together.
/// @param m1 first matrix.
/// @param m2 second matrix.
/// @return The sum of the matrices.
vector<vector<double>> sum_matrix(vector<vector<double>> m1,
                                  const vector<vector<double>> m2)
{
    //   NB: in the future it may be important to use the compiler defined sizes
    //   for platform portability
    size_t d1 = static_cast<size_t>(m1.size()),
           d2 = static_cast<size_t>(m1[0].size());

    // If the matrices do not have the same dimensions we return a dimension error
    // for the GUI to handle
    if (d1 != static_cast<size_t>(m2.size()) || d2 != static_cast<size_t>(m2[0].size()))
        return d_err;

    for (size_t i = 0; i < d1; ++i)
        for (size_t j = 0; j < d2; ++j)
            m1[i][j] += m2[i][j];
    return m1;
}

vector<vector<double>> sub_matrix(vector<vector<double>> m1,
                                  const vector<vector<double>> m2)
{
    //   NB: in the future it may be important to use the compiler defined sizes
    //   for platform portability
    const int d1 = static_cast<int>(m1.size()),
              d2 = static_cast<int>(m1[0].size());

    // If the matrices do not have the same dimensions we return a dimension error
    // for the GUI to handle
    if (d1 != static_cast<int>(m2.size()) || d2 != static_cast<int>(m2[0].size()))
        return d_err;

    for (int i = 0; i < d1; ++i)
        for (int j = 0; j < d2; ++j)
            m1[i][j] -= m2[i][j];

    return m1;
}

/// @brief Creates a random matrix of given dimensions filled with values in
/// specified range.
/// @param d1 Number of rows.
/// @param d2 Number of columns.
/// @param min Lower bound of values.
/// @param max Upper bound of values.
/// @return The generated matrix.
vector<vector<double>> generate_random_matrix(const int d1, const int d2,
                                              const double min,
                                              const double max)
{
    // must be a 1x1 matrix at the minimum
    if (d1 < 1 || d2 < 1)
        return d_err;

    vector<vector<double>> matrix;
    vector<double> line;

    random_device rd;
    default_random_engine eng(rd());
    uniform_real_distribution<double> distr(min, max);

    for (int i = 0; i < d1; ++i)
    {
        for (int j = 0; j < d2; ++j)
            line.push_back(distr(eng));
        matrix.push_back(line);
        line.clear();
    }

    return matrix;
}

/// @brief Multiplies two matrices together.
/// @param m1 First matrix.
/// @param m2 Second matrix.
/// @return Product matrix.
vector<vector<double>> mult_matrix(const vector<vector<double>> m1,
                                   const vector<vector<double>> m2)
{
    const int r1 = static_cast<int>(m1.size()),
              c1 = static_cast<int>(m1[0].size()),
              r2 = static_cast<int>(m2.size()),
              c2 = static_cast<int>(m2[0].size());
    //   columns of first matrix must equal rows of second
    if (c1 != r2)
        return d_err;

    vector<vector<double>> m3(r1);

    // Try using the std::par_unseq execution policy from C++20 (it's actually
    // pretty cool how it implements SIMD)
    for (auto i = 0; i < r1; i++)
        m3[i] = vector<double>(c2, 0);

    for (int i = 0; i < r1; ++i)
        for (int j = 0; j < c2; ++j)
            for (int k = 0; k < r2; ++k)
                m3[i][j] += m1[i][k] * m2[k][j];

    return m3;
}

/// @brief Multiplies a matrix by a scalar value.
/// @param matrix The matrix to be multiplied.
/// @param scalar The scalar value.
/// @return The resulting matrix after scalar multiplication.
vector<vector<double>> scalar_multiply(const vector<vector<double>> matrix, const double scalar)
{
    const int rows = static_cast<int>(matrix.size());
    const int cols = static_cast<int>(matrix[0].size());

    vector<vector<double>> result(rows, vector<double>(cols, 0.0));

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result[i][j] = matrix[i][j] * scalar;

    return result;
}


/// @brief Scales the matrix upwards by a given constant (i.e. multiply every
/// value in the matrix).
/// @param m1 The matrix.
/// @param s Scaling constant.
/// @return The updated matrix.
vector<vector<double>> scale_up(vector<vector<double>> m1, const double s)
{
    const int d1 = static_cast<int>(m1.size()),
              d2 = static_cast<int>(m1[0].size());
    for (auto i = 0; i < d1; ++i)
        for (auto j = 0; j < d2; ++j)
            m1[i][j] *= s;
    return m1;
}

/// @brief Scales the matrix downwards by a given constant (i.e. divides every
/// value in the matrix).
/// @param m1 The matrix.
/// @param s Scaling constant.
/// @return The update matrix.
vector<vector<double>> scale_down(vector<vector<double>> m1, const double s)
{
    // cannot divide by zero
    // AS 2023: changed it to s == 0 (instead of !s) to prevent potential float-
    // 	      int point precision errs
    if (s == 0.0)
        return a_err;

    const int d1 = static_cast<int>(m1.size()),
              d2 = static_cast<int>(m1[0].size());
    for (auto i = 0; i < d1; ++i)
        for (auto j = 0; j < d2; ++j)
            m1[i][j] /= s;
    return m1;
}

/// @brief Transposes the matrix.
/// @param m1 Input matrix.
/// @return Transposed matrix.
vector<vector<double>> transpose(const vector<vector<double>> m1)
{
    const int d1 = static_cast<int>(m1.size()),
              d2 = static_cast<int>(m1[0].size());

    vector<vector<double>> m2(d2);

    for (auto i = 0; i < d2; ++i)
        m2[i] = vector<double>(d1);

    for (auto i = 0; i < d1; ++i)
        for (auto j = 0; j < d2; ++j)
            m2[j][i] = m1[i][j];

    return m2;
}

/// @brief Saves the current working set of matrices.
/// @param matrices All of the matrices in the current process.
/// @param filename File to save matrices.
/// @return True on success, false otherwise.
bool save_file(const vector<vector<vector<double>>> matrices,
               const char *filename)
{
    fstream f;
    f.open(filename, fstream::out | fstream::trunc);
    if (f.fail())
    {
        cout << "error (" << errno << "): error opening file" << endl;
        return false;
    }
    auto size = matrices.size();
    f << size;
    f << '\n';

    size_t d1, d2;
    for (auto matrix : matrices)
    {
        d1 = matrix.size();
        d2 = matrix[0].size();
        f << d1;
        f << ',';
        f << d2;
        f << '\n';
        for (size_t j = 0; j < d1; ++j)
            for (size_t k = 0; k < d2; ++k)
            {
                f << matrix[j][k];
                if (j != d1 - 1 || k != d2 - 1)
                    f << ',';
                else if (j == d1 - 1 && k == d2 - 1 && size != 1)
                    f << '\n';
            }
        size--;
    }
    return true;
}

/// @brief Swaps the rows of a matrix.
/// @param m The input matrix.
/// @param i First row.
/// @param j Second row.
void swap_row(vector<vector<double>> &m, const int i, const int j)
{
    // NB: the matrix is passed by reference to avoid copy on call since matrix
    // could be large which could hurt performance.
    const int n = static_cast<int>(m.size());
    for (int k = 0; k <= n; k++)
    {
        const double temp = m[i][k];
        m[i][k] = m[j][k];
        m[j][k] = temp;
    }
}

/// @brief Performs forward elimination on the matrix.
/// @param m The input matrix.
/// @return -1 on failure, positive constant otherwise.
int forward_elimination(vector<vector<double>> &m)
{
    const int n = static_cast<int>(m.size());
    for (int k = 0; k < n; k++)
    {
        int i_max = k;
        int v_max = m[i_max][k];

        for (int i = k + 1; i < n; i++)
            if (abs(m[i][k]) > v_max)
                v_max = m[i][k], i_max = i;

        if (!m[k][i_max])
            return k;

        if (i_max != k)
            swap_row(m, k, i_max);

        for (int i = k + 1; i < n; i++)
        {
            double f = m[i][k] / m[k][k];

            for (int j = k + 1; j <= n; j++)
                m[i][j] -= m[k][j] * f;

            m[i][k] = 0;
        }
    }
    return -1;
}

/// @brief Performs a backwards substition for a given matrix.
/// @param m Input matrix.
/// @return the result.
vector<double> backward_substitution(vector<vector<double>> &m)
{
    const int n = static_cast<int>(m.size());
    vector<double> x(n);

    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = m[i][n];

        for (int j = i + 1; j < n; j++)
            x[i] -= m[i][j] * x[j];

        x[i] = x[i] / m[i][i];
    }

    vector<double> res(begin(x), end(x));
    return res;
}

// gaussian elimination with partial pivoting
// returns true if successful, false if A is singular
// Modifies both A and b to store the results
bool gaussian_elimination(std::vector<std::vector<double>> &A, std::vector<double> &b)
{
    const size_t n = A.size();

    for (size_t i = 0; i < n; i++)
    {

        // find pivot row and swap
        size_t max = i;
        for (size_t k = i + 1; k < n; k++)
            if (abs(A[k][i]) > abs(A[max][i]))
                max = k;
        std::swap(A[i], A[max]);
        std::swap(b[i], b[max]);

        // singular or nearly singular
        if (abs(A[i][i]) <= 1e-10)
            return false;

        // pivot within A and b
        for (size_t k = i + 1; k < n; k++)
        {
            double t = A[k][i] / A[i][i];
            for (size_t j = i; j < n; j++)
            {
                A[k][j] -= A[i][j] * t;
            }
            b[k] -= b[i] * t;
        }
    }

    // back substitution
    for (int i = n - 1; i >= 0; i--)
    {
        for (size_t j = i + 1; j < n; j++)
            b[i] -= A[i][j] * b[j];
        b[i] = b[i] / A[i][i];
    }
    return true;
}

// Perform QR factorization of A into Q and R matrices
// A is an m x n matrix with m >= n
// Q is an m x n orthogonal matrix
// R is an n x n upper-triangular matrix
pair<vector<vector<double>>, vector<vector<double>>> qr_factorization(vector<vector<double>> A)
{
    const int m = A.size();
    const int n = A[0].size();

    // Initialize Q and R matrices with zeros
    vector<vector<double>> Q(m, vector<double>(n, 0.0));
    vector<vector<double>> R(n, vector<double>(n, 0.0));

    // Calculate R matrix using Gram-Schmidt orthogonalization
    for (int j = 0; j < n; j++)
    {
        // Orthogonalize the jth column of A against previous columns
        for (int i = 0; i < m; i++)
        {
            Q[i][j] = A[i][j];
            for (int k = 0; k < j; k++)
            {
                double dot = 0;
                for (int l = 0; l < m; l++)
                {
                    dot += Q[l][k] * A[l][j];
                }
                Q[i][j] -= Q[i][k] * dot;
            }
        }
        // Calculate the jth row of R
        for (int i = 0; i < j; i++)
        {
            double dot = 0;
            for (int k = 0; k < m; k++)
            {
                dot += Q[k][i] * A[k][j];
            }
            R[i][j] = dot;
        }
        double norm = 0;
        for (int i = 0; i < m; i++)
        {
            norm += Q[i][j] * Q[i][j];
        }
        norm = sqrt(norm);
        for (int i = 0; i < m; i++)
        {
            Q[i][j] /= norm;
        }
        R[j][j] = norm;
    }
    return make_pair(Q, R);
}

// Perform LU factorization with partial pivoting on the input matrix A
// The factorization is stored in place in A as lower and upper triangular matrices
// with the diagonal of L stored as ones
// Return the permutation matrix as a vector of indices
std::vector<int> lu_factorization_inplace(std::vector<std::vector<double>> &A)
{
    const size_t n = A.size();
    std::vector<int> p(n);
    if (n != A[0].size())
    {
        throw invalid_argument("Error: Matrix must be square nxn");
    }
    // Initialize the permutation matrix to the identity matrix
    for (size_t i = 0; i < n; i++)
    {
        p[i] = i;
    }

    // Perform LU factorization with partial pivoting
    for (size_t k = 0; k < n - 1; k++)
    {
        // Find pivot row and swap
        size_t pivot_row = k;
        double pivot_val = std::abs(A[k][k]);
        for (size_t i = k + 1; i < n; i++)
        {
            double val = std::abs(A[i][k]);
            if (val > pivot_val)
            {
                pivot_val = val;
                pivot_row = i;
            }
        }
        if (pivot_val == 0)
        {
            throw std::runtime_error("Singular matrix");
        }
        if (pivot_row != k)
        {
            std::swap(p[k], p[pivot_row]);
            std::swap(A[k], A[pivot_row]);
        }

        // Eliminate entries below the pivot
        for (size_t i = k + 1; i < n; i++)
        {
            double factor = A[i][k] / A[k][k];
            A[i][k] = factor;
            for (size_t j = k + 1; j < n; j++)
            {
                A[i][j] -= factor * A[k][j];
            }
        }
    }

    return p;
}

tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>> lu_factorization(std::vector<std::vector<double>> &A)
{
    vector<vector<double>> L(A.size(), vector<double>(A.size(), 0.0));
    vector<vector<double>> U(A.size(), vector<double>(A.size(), 0.0));
    vector<vector<double>> P(A.size(), vector<double>(A.size(), 0.0));
    for (size_t i = 0; i < A.size(); i++)
    {
        for (size_t j = 0; j < A.size(); j++)
        {
            L[i][j] = A[i][j];
        }
    }
    vector<int> p = lu_factorization_inplace(L);
    for (size_t i = 0; i < A.size(); i++)
    {
        for (size_t j = 0; j < A.size(); j++)
        {
            if (i == j)
            {
                U[i][j] = L[i][j];
                L[i][j] = 1;
            }
            else if (i < j)
            {
                U[i][j] = L[i][j];
                L[i][j] = 0;
            }
            else
            {
                U[i][j] = 0;
            }
        }
    }
    for (size_t i = 0; i < A.size(); i++)
    {
        P[i][p[i]] = 1;
    }
    return std::make_tuple(P, L, U);
}

// Perfome Cholesky factorization on the input matrix A
// Return the lower triangular matrix L

vector<vector<double>> cholesky_factorization(std::vector<std::vector<double>> &A)
{
    const size_t n = A.size();
    vector<vector<double>> L(n, vector<double>(n, 0.0));
    if (n != A[0].size())
    {
        throw invalid_argument("Error: Matrix must be square nxn");
    }
    // Perform Cholesky factorization
    for (size_t j = 0; j < n; j++)
    {
        double sum = 0.0;
        for (size_t k = 0; k < j; k++)
        {
            sum += L[j][k] * L[j][k];
        }
        double d = A[j][j] - sum;
        if (d < 0.0)
        {
            throw invalid_argument("Error: Matrix is not positive definite");
        }
        L[j][j] = sqrt(d);
        for (size_t i = j + 1; i < n; i++)
        {
            double sum = 0.0;
            for (size_t k = 0; k < j; k++)
            {
                sum += L[i][k] * L[j][k];
            }
            L[i][j] = (A[i][j] - sum) / L[j][j];
        }
    }
    return L;
}

// Perform LDL^T factorization on the input matrix A
pair<std::vector<std::vector<double>>, std::vector<double>> ldlt_factorization(std::vector<std::vector<double>> &A)
{
    const int n = A.size();

    // Initialize D to the diagonal of A and L to the identity matrix
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0));
    std::vector<double> D(n, 0);
    for (int i = 0; i < n; i++)
    {
        D[i] = A[i][i];
        L[i][i] = 1;
    }

    // Perform factorization
    for (int j = 0; j < n; j++)
    {
        double sum = 0;
        for (int i = 0; i < j; i++)
        {
            sum += L[j][i] * D[i] * L[j][i];
        }
        D[j] -= sum;

        for (int i = j + 1; i < n; i++)
        {
            sum = 0;
            for (int k = 0; k < j; k++)
            {
                sum += L[i][k] * D[k] * L[j][k];
            }
            L[i][j] = (A[i][j] - sum) / D[j];
        }
    }

    return make_pair(L, D);
}

// Solves Ax = b using Gauss-Seidel method
// A is the matrix and b is the right-hand side vector
// x is the initial guess for the solution and max_iter is the maximum number of iterations
// Returns true if successful, false otherwise
std::vector<double> gauss_seidel(const std::vector<std::vector<double>> &A, const std::vector<double> &b, const double tol,const int max_iter)
{
    const int n = A.size();

    std::vector<double> x(n, 0.0);
    std::vector<double> x_new(n, 0.0);

    int iter = 0;
    double diff = tol + 1.0;

    // Run Gauss-Seidel iterations
    while (iter < max_iter && diff > tol)
    {
        for (int i = 0; i < n; i++)
        {
            double sum1 = 0;
            double sum2 = 0;
            for (int j = 0; j < i; j++)
            {
                sum1 += A[i][j] * x_new[j];
            }
            for (int j = i + 1; j < n; j++)
            {
                sum2 += A[i][j] * x[j];
            }
            x_new[i] = (b[i] - sum1 - sum2) / A[i][i];
        }
        // Calculate the norm of the difference between x and x_new
        diff = 0.0;
        for (int i = 0; i < n; i++)
        {
            double abs_diff = std::abs(x_new[i] - x[i]);
            if (abs_diff > diff)
            {
                diff = abs_diff;
            }
        }
        x = x_new;
        iter++;
    }
    cerr << "Gauss_seidel: "<< iter <<endl;
    return x;
}

// Solve the linear system Ax = b using Jacobi iteration
// A is a square matrix of size n x n
// b is a vector of length n
// tol is the tolerance for convergence
// max_iter is the maximum number of iterations to perform
// Return the solution x as a vector of length n
std::vector<double> jacobi_iteration(const std::vector<std::vector<double>> &A,
                                     const std::vector<double> &b,
                                     const double tol,
                                     const int max_iter)
{
    const int n = A.size();
    std::vector<double> x(n, 0.0);
    std::vector<double> x_new(n, 0.0);

    int iter = 0;
    double diff = tol + 1.0;
    while (iter < max_iter && diff > tol)
    {
        // Compute new iterate
        for (int i = 0; i < n; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
            {
                if (j != i)
                {
                    sum += A[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - sum) / A[i][i];
        }
        // Calculate the norm of the difference between x and x_new
        diff = 0.0;
        for (int i = 0; i < n; i++)
        {
            double abs_diff = std::abs(x_new[i] - x[i]);
            if (abs_diff > diff)
            {
                diff = abs_diff;
            }
        }
        x = x_new;
        iter++;
    }
    cerr << "Jacobi Itertions: "<< iter <<endl;
    return x;
}


/**
 * @brief Solve the linear system Ax = b using SSOR iteration with relaxation parameter omega
 * A is a square matrix of size n x n
 * b is a vector of length n
 * tol is the tolerance for convergence
 * max_iter is the maximum number of iterations to perform
 * omega is the relaxation parameter (0 < omega < 2)
 * Return the solution x as a vector of length n
 * 
 * @param A 
 * @param b 
 * @param tol 
 * @param max_iter 
 * @param omega 
 * @return std::vector<double> 
 */
std::vector<double> ssor_iteration(const std::vector<std::vector<double>> &A,
                                   const std::vector<double> &b,
                                   const double tol,
                                   const int max_iter,
                                   const double omega)
{
    const int n = A.size();
    std::vector<double> x(n, 0.0);
    std::vector<double> x_new(n, 0.0);

    int iter = 0;
    double diff = tol + 1.0;

    while (iter < max_iter && diff > tol)
    {
        // Forward Sweep (Red-Black Ordering)
        for (int i = 0; i < n; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
            {
                if (j != i)
                {
                    sum += A[i][j] * x_new[j];
                }
            }
            x_new[i] = (1.0 - omega) * x[i] + (omega / A[i][i]) * (b[i] - sum);
        }

        // Backward Sweep (Black-Red Ordering)
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
            {
                if (j != i)
                {
                    sum += A[i][j] * x_new[j];
                }
            }
            x_new[i] = (1.0 - omega) * x[i] + (omega / A[i][i]) * (b[i] - sum);
        }

        // Calculate the norm of the difference between x and x_new
        diff = 0.0;
        for (int i = 0; i < n; i++)
        {
            double abs_diff = std::abs(x_new[i] - x[i]);
            if (abs_diff > diff)
            {
                diff = abs_diff;
            }
        }

        x = x_new;
        iter++;
    }
    cerr << "SSOR Itertions : "<< iter <<endl;
    return x;
}


/**
 * @brief Incomplete Cholesky Factorization
 * 
 * @param A 
 * @param tol 
 * @return std::vector<vector<double>> 
 */
std::vector<vector<double>> incompleteCholesky(const vector<vector<double>> &A, double tol)
{
    int n = A.size();
    vector<vector<double>> L(n, vector<double>(n));

    // Compute the lower triangle of A
    for (int i = 0; i < n; ++i)
    {
        // cout << "A[i][i]: " << A[i][i] << endl;
        double d = A[i][i];
        for (int k = 0; k < i; ++k)
        {
            // cout << "i, k, d: " << i << ", " << k << ", " << d << endl;
            d -= L[i][k] * L[i][k];
        }
        // cout << "d: " << d << endl;
        L[i][i] = std::sqrt(std::max(d, 0.0));
        for (int j = i + 1; j < n; ++j)
        {
            double s = 0.0;
            for (int k = 0; k < i; ++k)
            {
                s += L[i][k] * L[j][k];
                // cout << "j, k, i, s: " << j << ", " << k << ", " << i << ", " << s << endl;
            }
            L[j][i] = (A[j][i] - s) / L[i][i];
            // d -= L[j][k] * L[j][k];
        }
    }

    return L;
}

/**
 * @brief Finding Matrix determinant, once matrix size gets to over 1000x1000, this will not compute in a fast
 * enough time, use determinant using LU Decomposition.
 * 
 * Time complexity O(n!)
 * 
 * @param matrix 
 * @param numRows 
 * @param numCols 
 * @return double 
 */
double find_matrix_determinant(std::vector<std::vector<double>> &matrix) {
    int dimension = matrix.size();

    if (dimension == 1) {
        return matrix[0][0];
    }

    double determinant = 0;

    for (int i = 0; i < dimension; ++i) {
        std::vector<std::vector<double>> subMat;
        for (int j = 1; j < dimension; ++j) {
            std::vector<double> row;
            for (int k = 0; k < dimension; ++k) {
                if (k != i) {
                    row.push_back(matrix[j][k]);
                }
            }
            subMat.push_back(row);
        }
        determinant += (i % 2 == 0 ? 1 : -1) * matrix[0][i] * find_matrix_determinant(subMat);
    }

    return determinant;
}

/**
 * @brief Matrix determinant using LU decomposition. Using boiler plate LU inplace then finding
 * det(A) = det(L) * det(U).
 * 
 * det(U) and det(L) is the product along the diagonals
 * 
 * Time complexity O(n^3)
 * 
 * @param A 
 * @return double 
 */
double matrix_determinant_lu(std::vector<std::vector<double>> &A)
{
    int n = A.size();
    vector<vector<double>> L(A.size(), vector<double>(A.size(), 0.0));
    vector<vector<double>> U(A.size(), vector<double>(A.size(), 0.0));
    vector<vector<double>> P(A.size(), vector<double>(A.size(), 0.0));
    double determinant = -1.0;
    for (size_t i = 0; i < A.size(); i++)
    {
        for (size_t j = 0; j < A.size(); j++)
        {
            L[i][j] = A[i][j];
        }
    }
    vector<int> p = lu_factorization_inplace(L);
    for (size_t i = 0; i < A.size(); i++)
    {
        for (size_t j = 0; j < A.size(); j++)
        {
            if (i == j)
            {
                U[i][j] = L[i][j];
                L[i][j] = 1;
            }
            else if (i < j)
            {
                U[i][j] = L[i][j];
                L[i][j] = 0;
            }
            else
            {
                U[i][j] = 0;
            }
        }
    }
    for (size_t i = 0; i < A.size(); i++)
    {
        P[i][p[i]] = 1;
    }
    double u_determinant = 1.0;
    for (int i = 0; i < n; i++) {
        u_determinant *= U[i][i];
    }

    double l_determinant = 1.0;
    for (int i = 0; i < n; ++i) {
        l_determinant *= L[i][i];
    }  

    determinant = determinant * u_determinant * l_determinant;

    return determinant;
    
}

std::vector<std::vector<double>> construct_identity_matrix(int rows, int columns) {
    std::vector<std::vector<double>> matrix(rows, std::vector<double>(columns, 0.0));
     for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            if (i == j) {
                matrix[i][j] = 1.0; 
            } else {
                matrix[i][j] = 0.0; 
            }
        }
    }
    return matrix;
}

/**
 * @brief Matrix inverse using Guass Elimination (Diep) using an augmented version of A
 * 
 * @param A 
 * @return std::vector<std::vector<double>> 
 */
bool matrix_inverse(std::vector<std::vector<double>> &A) {
    /* TODO: add checks to ensure that it is invertible */
    /* temporary -- can cause n! */
    if(find_matrix_determinant(A) == 0) {
        return false; /* matrix is not invertable */
    }
    const size_t n = A.size();

    // Agugment A with the Identity Matrix so that we can do Guass and find the inverse
    std::vector<std::vector<double>> augmentedA(n, std::vector<double>(2 * n));
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            augmentedA[i][j] = A[i][j];
            if (i == j) {
                augmentedA[i][j+n] = 1.0;
            } else {
                augmentedA[i][j+n] = 0.0;
            }
        }
    }

    for (size_t i = 0; i < n; i++) {
        size_t max = i;
        for (size_t k = i + 1; k < n; k++) {
            if (abs(augmentedA[k][i]) > abs(augmentedA[max][i])) {
                max = k;
            }
        }
        std::swap(augmentedA[i], augmentedA[max]);
        for (size_t k = 0; k < n; k++) {
            if (k != i) {
                double factor = augmentedA[k][i] / augmentedA[i][i];
                for (size_t j = 0; j < 2 * n; j++) {
                    augmentedA[k][j] -= factor * augmentedA[i][j];
                }
            }
        }

        double pivot = augmentedA[i][i];
        for (size_t j = 0; j < 2 * n; j++) {
            if(pivot != 0.0){
                augmentedA[i][j] /= pivot;
            }else{  
                return false;
            }
        }
    }

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            A[i][j] = augmentedA[i][j + n];
        }
    }
    return true; /* result is stored in A */
}

// dot product of two vectors
double dot(const std::vector<double> &x, const std::vector<double> &y)
{
    const size_t n = x.size();
    double sum = 0.0;
    for (size_t i = 0; i < n; i++)
    {
        sum += x[i] * y[i];
    }
    return sum;
}

// left multiply a matrix by a vector
std::vector<double> left_mult_vector(const std::vector<std::vector<double>> &A, const std::vector<double> &x)
{
    const size_t n = A.size();
    std::vector<double> y(n, 0.0);
    for (size_t i = 0; i < n; i++)
    {
        y[i] = dot(A[i], x);
    }
    return y;
}

// GCR method solves Ax = b
// A is the matrix and b is the right-hand side vector
// x is the initial guess for the solution
// max_iter is the maximum number of iterations to perform
// tol is the tolerance for convergence
// Return the solution x as a vector of length n
std::vector<double> gcr(const std::vector<std::vector<double>> &A,
                        const std::vector<double> &b,
                        const double tol,
                        const int max_iter)

{

    const size_t n = A.size();
    std::vector<double> x(n, 0.0);
    std::vector<double> r(n, 0);
    std::vector<std::vector<double>> p, beta;
    // r = b - Ax
    r = b;
    for (size_t i = 0; i < n; i++)
    {
        r[i] -= dot(A[i], x);
    }
    // p_0 = r
    p.push_back(r);
    for (int j = 0; j < max_iter; j++)
    {   
        // Check for convergence
        double norm_r = sqrt(dot(r, r));
        if (norm_r < tol)
        {
            cerr << "GCR Itertions : "<< j <<endl;
            return x;
        }
        // A*p
        std::vector<double> Ap = left_mult_vector(A, p[j]);
        // alpha = (r, Ap) / (Ap, Ap)
        double alpha = dot(r, Ap) / dot(Ap, Ap);
        // x = x + alpha*p and r = r - alpha*Ap
        // std::vector<double> x_new = x;
        for (size_t i = 0; i < n; i++)
        {
            x[i] += alpha * p[j][i];
            r[i] -= alpha * Ap[i];
        }
        // Compute beta_ij
        std::vector<double> beta;
        for (int i=0; i < j; i++) {
            std::vector<double> Ar = left_mult_vector(A, r);
            std::vector<double> Api = left_mult_vector(A, p[i]);
            // beta_ij = (r_new, A*p_i) / (A*p_i, A*p_i)
            beta.push_back(dot(Ar, Api) / dot(Api, Api));
        }

        std::vector<double> p_new = r;
        for (int i=0; i < j; i++) {
            for (size_t k=0; k < n; k++) {
                p_new[k] -= beta[i] * p[i][k];
            }
        }
        p.push_back(p_new);
    }
    cerr << "GCR Itertions : "<< max_iter <<endl;
    return x;
}


/// @brief Creates dense matrix for .mtx file
/// @tparam T the type of matrix
/// @param fileName the name of the file to import
/// @return a new dense matrix from the filename
template <typename T>
std::vector<std::vector<T>> load_fileMatrix(const std::string& fileName) {
    std::ifstream file(fileName);

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << fileName << std::endl;
        exit(1); // Handle the error as per your requirements
    }

    int num_row = 0, num_col = 0, num_lines = 0;

    // Ignore comments and headers starting with '%'
    while (file.peek() == '%')
        file.ignore(2048, '\n');

    // Read number of rows, columns, and non-zero values
    file >> num_row >> num_col >> num_lines;

    std::vector<std::vector<T>> matrix(num_row, std::vector<T>(num_col, 0));

    T data;
    int row, col;

    for (int i = 0; i < num_lines; i++) {
        file >> row >> col >> data;
        row--; // Convert from 1-based to 0-based index
        col--; // Convert from 1-based to 0-based index

        matrix[row][col] = data;
    }

    file.close();

    return matrix;
}

#include <iostream>
#include <vector>
#include <cassert>

#include <chrono>
class timer {
public:
    std::chrono::time_point<std::chrono::high_resolution_clock> lastTime;
    timer() : lastTime(std::chrono::high_resolution_clock::now()) {}
    inline double elapsed() {
        std::chrono::time_point<std::chrono::high_resolution_clock> thisTime=std::chrono::high_resolution_clock::now();
        double deltaTime = std::chrono::duration<double>(thisTime-lastTime).count();
        lastTime = thisTime;
        return deltaTime;
    }
};
#include <iomanip>

void printMatrix(vector<vector<double>> &A) {
    // Set the width of each output element to 8 characters
    const int width = 2;
    // Loop over each row of the matrix
    for (size_t i = 0; i < A.size(); i++) {
        // Loop over each element in the row
        for (size_t j = 0; j < A[i].size(); j++) {
            // Set the output width for the element
            cout << setw(width) << A[i][j];
        }
        // End the row with a newline
        cout << endl;
    }
}

void printVector(vector<double> &B) {
    // Set the width of each output element to 8 characters
    const int width = 2;
        // Loop over each element in the row
        for (size_t j = 0; j < B.size(); j++) {
            // Set the output width for the element
            cout << setw(width) << B[j];
        }
        // End the row with a newline
        cout << endl;
}

// int main() {
//     // size_t N = 10;
//     // vector<vector<double>> A(N, vector<double>(N, 0.0));
//     // vector<double> B(N, 0.0);
//     // B[0] = 1.0;
//     // for (int i = 0; i < N; ++i) {
//     //     A[i][i] = 2.0;
//     //     if (i >0) {
//     //         A[i][i-1] = -1.0;
//     //     } 
//     //     if (i < N - 1) {
//     //         A[i][i+1] = -1.0;
//     //     }
//     // }
//     // printMatrix(A);

//     // cout << endl;
//     // printVector(B);
//     // vector<double> jacobi_solved = jacobi_iteration(A,B,1e-6,1000);
//     // printVector(jacobi_solved);
//     // vector<double> ssor_solved = ssor_iteration(A,B,1e-6,1000,0.5);
//     // printVector(ssor_solved);
//     // gaussian_elimination(A,B);
//     // printVector(B);

//     timer stopwatch;
//     std::vector<double> jacobi;
//     std::vector<double> omega_5;
//     std::vector<double> omega_10;
//     std::vector<double> gauss_seidel_time;
//     std::vector<double> gauss_time_vector;
//     std::vector<int> size = {500,1000,1500,2000,2500,3000,3500,4000};
//     //std::vector<int> size = {500,600,700};
//     for(size_t i = 0 ; i <size.size();i++){
//         size_t N = size[i];
//         vector<vector<double>> A(N, vector<double>(N, 0.0));
//         vector<double> B(N, 0.0);
//         B[0] = 1.0;
//         for (int i = 0; i < N; ++i) {
//             A[i][i] = 2.0;
//             if (i >0) {
//                 A[i][i-1] = -1.0;
//             } 
//             if (i < N - 1) {
//                 A[i][i+1] = -1.0;
//             }
//         }
//         stopwatch.elapsed();
//         jacobi_iteration(A,B,1e-6,1000);
//         jacobi.push_back(stopwatch.elapsed());
//         stopwatch.elapsed();
//         ssor_iteration(A,B,1e-6,1000,1.0);
//         omega_10.push_back(stopwatch.elapsed());
//         stopwatch.elapsed();
//         gauss_seidel(A,B,1e-6,1000);
//         gauss_seidel_time.push_back(stopwatch.elapsed());
//         stopwatch.elapsed();
//         gaussian_elimination(A,B);
//         gauss_time_vector.push_back(stopwatch.elapsed());
//     }
//     cerr<< "Jacobi: ";
//     for(auto val : jacobi){
//         cerr<< val << ",";
//     }
//     cerr<< endl;
//     cerr<< "SSOR 1.0: ";
//     for(auto val : omega_10){
//         cerr<< val << ",";
//     }
//     cerr<< endl;
//     cerr<< "Gauss Seidel: ";
//     for(auto val : gauss_seidel_time){
//         cerr<< val << ",";
//     }
//     cerr<< endl;
//     cerr<< "Gauss Elimination: ";
//     for(auto val : gauss_time_vector){
//         cerr<< val << ",";
//     }
//     cerr<< endl;

//   return 0;
// }

// int main() {
//     std::string to_load = "bcsstk10.mtx";
//     vector<vector<double>> A = load_fileMatrix<double>(to_load);
//     vector<double> B(A.size(), 0);
//     B[0] = 1.0;
//     timer stopwatch;
//     stopwatch.elapsed();
//     jacobi_iteration(A,B,1e-16,1000);
//     cerr << "Jacobi: "<< stopwatch.elapsed() <<endl;
//     stopwatch.elapsed();
//     gauss_seidel(A,B,1e-16,1000);
//     cerr << "Gauss_seidel: "<< stopwatch.elapsed() <<endl;
//     stopwatch.elapsed();
//     ssor_iteration(A,B,1e-16,1000,0.5);
//     cerr << "SSOR w = 0.5: "<< stopwatch.elapsed() <<endl;
//     stopwatch.elapsed();
//     gaussian_elimination(A,B);
//     cerr << "Gauss Elimination: "<< stopwatch.elapsed() <<endl;

//     // vector<double> jacobi_solved = jacobi_iteration(A,B,1e-6,1000);
//     //printVector(jacobi_solved);
//     // cerr <<endl;
//     // gaussian_elimination(A,B);
//     //printVector(B);
//     // cerr <<endl;
//     return 0;
// }

//     // cout << endl;
//     // printVector(B);
//     // vector<double> jacobi_solved = jacobi_iteration(A,B,1e-6,1000);
//     // printVector(jacobi_solved);
//     // vector<double> ssor_solved = ssor_iteration(A,B,1e-6,1000,0.5);
//     // printVector(ssor_solved);
//     // gaussian_elimination(A,B);
//     // printVector(B);

//     timer stopwatch;
//     std::vector<double> jacobi;
//     std::vector<double> omega_5;
//     std::vector<double> omega_10;
//     std::vector<double> gauss_seidel_time;
//     std::vector<double> gauss_time_vector;
//     std::vector<int> size = {500,1000,1500,2000,2500,3000,3500,4000};
//     //std::vector<int> size = {500,600,700};
//     for(size_t i = 0 ; i <size.size();i++){
//         size_t N = size[i];
//         vector<vector<double>> A(N, vector<double>(N, 0.0));
//         vector<double> B(N, 0.0);
//         B[0] = 1.0;
//         for (int i = 0; i < N; ++i) {
//             A[i][i] = 2.0;
//             if (i >0) {
//                 A[i][i-1] = -1.0;
//             } 
//             if (i < N - 1) {
//                 A[i][i+1] = -1.0;
//             }
//         }
//         stopwatch.elapsed();
//         jacobi_iteration(A,B,1e-6,1000);
//         jacobi.push_back(stopwatch.elapsed());
//         stopwatch.elapsed();
//         ssor_iteration(A,B,1e-6,1000,1.0);
//         omega_10.push_back(stopwatch.elapsed());
//         stopwatch.elapsed();
//         gauss_seidel(A,B,1e-6,1000);
//         gauss_seidel_time.push_back(stopwatch.elapsed());
//         stopwatch.elapsed();
//         gaussian_elimination(A,B);
//         gauss_time_vector.push_back(stopwatch.elapsed());
//     }
//     cerr<< "Jacobi: ";
//     for(auto val : jacobi){
//         cerr<< val << ",";
//     }
//     cerr<< endl;
//     cerr<< "SSOR 1.0: ";
//     for(auto val : omega_10){
//         cerr<< val << ",";
//     }
//     cerr<< endl;
//     cerr<< "Gauss Seidel: ";
//     for(auto val : gauss_seidel_time){
//         cerr<< val << ",";
//     }
//     cerr<< endl;
//     cerr<< "Gauss Elimination: ";
//     for(auto val : gauss_time_vector){
//         cerr<< val << ",";
//     }
//     cerr<< endl;

//   return 0;
// }

// int main() {
//     // size_t N = 10;
//     // vector<vector<double>> A(N, vector<double>(N, 0.0));
//     // vector<double> B(N, 0.0);
//     // B[0] = 1.0;
//     // for (int i = 0; i < N; ++i) {
//     //     A[i][i] = 2.0;
//     //     if (i >0) {
//     //         A[i][i-1] = -1.0;
//     //     } 
//     //     if (i < N - 1) {
//     //         A[i][i+1] = -1.0;
//     //     }
//     // }
//     // printMatrix(A);

//     // cout << endl;
//     // printVector(B);
//     // vector<double> jacobi_solved = jacobi_iteration(A,B,1e-6,1000);
//     // printVector(jacobi_solved);
//     // vector<double> ssor_solved = ssor_iteration(A,B,1e-6,1000,0.5);
//     // printVector(ssor_solved);
//     // gaussian_elimination(A,B);
//     // printVector(B);

//     timer stopwatch;
//     std::vector<double> jacobi;
//     std::vector<double> omega_5;
//     std::vector<double> omega_10;
//     std::vector<double> gauss_seidel_time;
//     std::vector<double> gauss_time_vector;
//     std::vector<int> size = {500};
//     for(size_t i = 0 ; i <size.size();i++){
//         size_t N = size[i];
//         vector<vector<double>> A(N, vector<double>(N, 0.0));
//         vector<double> B(N, 0.0);
//         B[0] = 1.0;
//         for (int i = 0; i < N; ++i) {
//             A[i][i] = 2.0;
//             if (i >0) {
//                 A[i][i-1] = -1.0;
//             } 
//             if (i < N - 1) {
//                 A[i][i+1] = -1.0;
//             }
//         }
//         jacobi = jacobi_iteration(A,B,1e-6,1000);
        
//     }
//     cerr<< "Jacobi: ";
//     for(auto val : jacobi){
//         cerr<< val << ",";
//     }
//     cerr<< endl;
//   return 0;
// }

// int main() {
//     timer stopwatch;
//     std::vector<double> dense;
//     std::vector<double> sparse;
//     // std::vector<int> size = {1000,1500,2000,2500,3000,3500,4000,4500};
//     std::vector<int> size = {500,600,700,800,900,1000};
//     for(size_t i = 0 ; i <size.size();i++){
//         std::vector<std::vector<double> > A = generate_random_matrix_sparse(size[i],size[i],1,10000,3);
//         CSRMatrix<double> B = from_vector(A);
//         stopwatch.elapsed();
//         mult_matrix(A,A);
//         dense.push_back(stopwatch.elapsed());
//         stopwatch.elapsed();
//         multiply_matrixCSR(B, B);
//         sparse.push_back(stopwatch.elapsed());
//     }
//     for(auto val : dense){
//         cerr<< val << ",";
//     }
//     cerr<< endl;
//     for(auto val : sparse){
//         cerr<< val << ",";
//     }
//     cerr<< endl;
//     return 0;
// }