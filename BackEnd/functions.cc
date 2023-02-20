// functions.cc
// Contains all of the functions for performing operations on matrices.

// Potential improvement, implement the functions using iterators instead of
// counter based loops.
// See: algorithm.h by GCC

// AS 12 February 2023
// 	- I removed the Open MP configurationse since they were mostly
// 	  redundant (e.g. local vars are by default private).
//	- Looking back it probably would have been better to have used Intel's 
//        One API's TBB (it is what it is).
//
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace std;

/// error with dimensions in a matrix
vector<vector<double> > d_err{{INT_MAX, INT_MAX, INT_MIN},
                              {INT_MIN, INT_MIN, INT_MAX}};
/// arithmetic error
vector<vector<double> > a_err{{INT_MIN, INT_MIN, INT_MAX},
                              {INT_MAX, INT_MAX, INT_MIN}};

/// @brief Reads a file and transforms the bytes into a 2d vector (a.k.a
/// matrix).
/// @param filename The name of the file with the data.
/// @return 2d vector containing the data.
vector<vector<vector<double> > > read_file(char *filename) {
  try {
    fstream file(filename);
    if (file.fail()) {
      cout << "error (" << errno << "): failed to open file '" << filename
           << "'" << endl;
      return vector<vector<vector<double> > >();
    }

    string row;
    vector<string> data;
    while (getline(file, row))
      data.push_back(row);

    string s_temp;
    vector<string> f_data;

    for (auto i : data) {
      stringstream stream(i);
      while (getline(stream, s_temp, ','))
        f_data.push_back(s_temp);
    }

    int num;

    int d1, d2;

    vector<vector<double> > temp;
    vector<vector<vector<double> > > p_matrices;

    data = f_data;

    vector<double> line;

    num = (int)stof(data[0]);
    data.erase(data.begin(), data.begin() + 1);

#pragma omp parallel
    for (int i = 0; i < num; ++i) {
      d1 = (int)stof(data[0]);
      d2 = (int)stof(data[1]);
      data.erase(data.begin(), data.begin() + 2);

      for (int j = 0; j < d1; ++j) {
        for (int k = 0; k < d2; ++k) {
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
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return {{}};
  }
}

/// @brief Adds to matrices together.
/// @param m1 first matrix.
/// @param m2 second matrix.
/// @return The sum of the matrices.
vector<vector<double> > sum_matrix(vector<vector<double> > m1,
                                   const vector<vector<double> > m2) {
  //   NB: in the future it may be important to use the compiler defined sizes
  //   for platform portability
  const int d1 = static_cast<int>(m1.size()),
            d2 = static_cast<int>(m1[0].size());

  // If the matrices do not have the same dimensions we return a dimension error
  // for the GUI to handle
  if (d1 != static_cast<int>(m2.size()) || d2 != static_cast<int>(m2[0].size()))
    return d_err;

#pragma omp parallel
  for (size_t i = 0; i < d1; ++i)
    for (size_t j = 0; j < d2; ++j)
      m1[i][j] -= m2[i][j];

  return m1;
}

vector<vector<double> > sub_matrix(vector<vector<double> > m1,
                                   const vector<vector<double> > m2) {
  //   NB: in the future it may be important to use the compiler defined sizes
  //   for platform portability
  const int d1 = static_cast<int>(m1.size()),
            d2 = static_cast<int>(m1[0].size());

  // If the matrices do not have the same dimensions we return a dimension error
  // for the GUI to handle
  if (d1 != static_cast<int>(m2.size()) || d2 != static_cast<int>(m2[0].size()))
    return d_err;

#pragma omp parallel
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
vector<vector<double> > generate_random_matrix(const int d1, const int d2,
                                               const double min,
                                               const double max) {
  // must be a 1x1 matrix at the minimum
  if (d1 < 1 || d2 < 1)
    return d_err;

  vector<vector<double> > matrix;
  vector<double> line;

  random_device rd;
  default_random_engine eng(rd());
  uniform_real_distribution<double> distr(min, max);

#pragma omp parallel
  for (int i = 0; i < d1; ++i) {
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
vector<vector<double> > mult_matrix(const vector<vector<double> > m1,
                                    const vector<vector<double> > m2) {
  const int r1 = static_cast<int>(m1.size()),
            c1 = static_cast<int>(m1[0].size()),
            r2 = static_cast<int>(m2.size()),
            c2 = static_cast<int>(m2[0].size());
  //   columns of first matrix must equal rows of second
  if (c1 != r2)
    return d_err;

  vector<vector<double> > m3(r1);

  // Try using the std::par_unseq execution policy from C++20 (it's actually
  // pretty cool how it implements SIMD)
#pragma omp parallel
  for (auto i = 0; i < r1; i++)
    m3[i] = vector<double>(c2, 0);

#pragma omp parallel
  for (int i = 0; i < r1; ++i)
    for (int j = 0; j < c2; ++j)
      for (int k = 0; k < r2; ++k)
        m3[i][j] += m1[i][k] * m2[k][j];

  return m3;
}

/// @brief Scales the matrix upwards by a given constant (i.e. multiply every
/// value in the matrix).
/// @param m1 The matrix.
/// @param s Scaling constant.
/// @return The updated matrix.
vector<vector<double> > scale_up(vector<vector<double> > m1, const double s) {
  const int d1 = static_cast<int>(m1.size()),
            d2 = static_cast<int>(m1[0].size());
#pragma omp parallel
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
vector<vector<double> > scale_down(vector<vector<double> > m1, const double s) {
  // cannot divide by zero
  // AS 2023: changed it to s == 0 (instead of !s) to prevent potential float-
  // 	      int point precision errs
  if (s == 0.0)
    return a_err;

  const int d1 = static_cast<int>(m1.size()),
            d2 = static_cast<int>(m1[0].size());
#pragma omp parallel
  for (auto i = 0; i < d1; ++i)
    for (auto j = 0; j < d2; ++j)
      m1[i][j] /= s;
  return m1;
}

/// @brief Transposes the matrix.
/// @param m1 Input matrix.
/// @return Transposed matrix.
vector<vector<double> > transpose(const vector<vector<double> > m1) {
  const int d1 = static_cast<int>(m1.size()),
            d2 = static_cast<int>(m1[0].size());

  vector<vector<double> > m2(d2);

#pragma omp parallel
  for (auto i = 0; i < d2; ++i)
    m2[i] = vector<double>(d1);

#pragma omp parallel
  for (auto i = 0; i < d1; ++i)
    for (auto j = 0; j < d2; ++j)
      m2[j][i] = m1[i][j];

  return m2;
}

/// @brief Saves the current working set of matrices.
/// @param matrices All of the matrices in the current process.
/// @param filename File to save matrices.
/// @return True on success, false otherwise.
bool save_file(const vector<vector<vector<double> > > matrices,
               const char *filename) {
  fstream f;
  f.open(filename, fstream::out | fstream::trunc);
  if (f.fail()) {
    cout << "error (" << errno << "): error opening file" << endl;
    return false;
  }
  auto size = matrices.size();
  f << size;
  f << '\n';

  size_t d1, d2;
#pragma omp parallel
  for (auto matrix : matrices) {
    d1 = matrix.size();
    d2 = matrix[0].size();
    f << d1;
    f << ',';
    f << d2;
    f << '\n';
    for (auto j = 0; j < d1; ++j)
      for (auto k = 0; k < d2; ++k) {
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
void swap_row(vector<vector<double> > &m, const int i, const int j) {
  // NB: the matrix is passed by reference to avoid copy on call since matrix
  // could be large which could hurt performance.
  const int n = static_cast<int>(m.size());
  for (int k = 0; k <= n; k++) {
    const double temp = m[i][k];
    m[i][k] = m[j][k];
    m[j][k] = temp;
  }
}

/// @brief Performs forward elimination on the matrix.
/// @param m The input matrix.
/// @return -1 on failure, positive constant otherwise.
int forward_elimination(vector<vector<double> > &m) {
  const int n = static_cast<int>(m.size());
  for (int k = 0; k < n; k++) {
    int i_max = k;
    int v_max = m[i_max][k];

    for (int i = k + 1; i < n; i++)
      if (abs(m[i][k]) > v_max)
        v_max = m[i][k], i_max = i;

    if (!m[k][i_max])
      return k;

    if (i_max != k)
      swap_row(m, k, i_max);

    for (int i = k + 1; i < n; i++) {
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
vector<double> backward_substitution(vector<vector<double> > &m) {
  const int n = static_cast<int>(m.size());
  vector<double> x(n);

  for (int i = n - 1; i >= 0; i--) {
    x[i] = m[i][n];

    for (int j = i + 1; j < n; j++)
      x[i] -= m[i][j] * x[j];

    x[i] = x[i] / m[i][i];
  }

  vector<double> res(begin(x), end(x));
  return res;
}

/// @brief Performs gaussian elmination.
/// @param m Input matrix
/// @return The computation result.
vector<vector<double> > gaussian_elimination(vector<vector<double> > &m) {
  const int singular_flag = forward_elimination(m);

  if (singular_flag != -1) {
    int n = static_cast<int>(m.size());
    if (m[singular_flag][n]) {
      cout << "inconsistent system"
           << "endl";
      vector<double> r(n + 1);
      return {r};
    }
    vector<double> r(n + 2);
    return {r};
  }
  // NB: the result needs to be wrapped in another vector since the program only
  // deals with matrices.
  return {backward_substitution(m)};
}

namespace py = pybind11;

PYBIND11_MODULE(functions, m) {
  m.doc() = "pybind11 example plugin";
  m.def("read_file", &read_file, "Read the contents of a file");
  m.def("sum_matrix", &sum_matrix, "Sum two matrices");
  m.def("sub_matrix", &sub_matrix, "Subtract two matrices");
  m.def("mult_matrix", &mult_matrix, "Multiply two matrices");
  m.def("scale_up", &scale_up, "Multiply a matrix by a constant factor");
  m.def("scale_down", &scale_down, "Divide a matrix by a constant factor");
  m.def("transpose", &transpose, "Tranpose a matrix");
  m.def("generate_random_matrix", &generate_random_matrix,
        "Generate a random matrix");
  m.def("save_file", &save_file, "save matrices to file");
  m.def("gaussian_elimination", &gaussian_elimination, "solve linear system");
}
