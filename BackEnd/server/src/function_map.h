#include <string>
#include <vector>

/* dense matrix operations */
std::vector<std::vector<double>> sum_matrix(const std::vector<std::vector<double>> m1, const std::vector<std::vector<double>> m2);
std::vector<std::vector<double>> scalar_multiply(const std::vector<std::vector<double>> matrix, const double scalar);
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> m1);
bool matrix_inverse(std::vector<std::vector<double>> &A);
bool gaussian_elimination(std::vector<std::vector<double>> &A, std::vector<double> &b);
std::vector<int> lu_factorization_inplace(std::vector<std::vector<double>> &A);
std::vector<double> jacobi_iteration(const std::vector<std::vector<double>> &A, const std::vector<double> &b, const double tol, const int max_iter);
bool gauss_seidel(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x, const int max_iter);

/* sparse matrix operations */