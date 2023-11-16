#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include "functions.cc"

using namespace std;

template <typename T>
class CSRMatrix
{
    // struct colValStruct {
    //     int col;
    //     T value;
    // };
public:
    size_t numRows, numColumns;
    // vector<colValStruct*> rowPointersVector;
    // fist value is col_ind and the second value is the value
    // vector<colValStruct> columnValueVector;
    vector<T> val;
    vector<size_t> col_ind;
    vector<size_t> row_ptr;
};
// TODO loadfile
// TODO savefile

/// @brief Gets a value from the compressed sparse row(CSR) matrix
/// @exception The search row and col must be less than dimensions of m1
/// @tparam T The type of the matrix
/// @param m1 The CSR matrix to get the value from
/// @param row The row of the value to get
/// @param col The column of the value to get
/// @return That value stored at row,column
template <typename T>
T get_matrixCSR(CSRMatrix<T> m1, size_t row, size_t col)
{
    if (m1.numRows <= row)
    {
        throw std::invalid_argument("The row being searched for is greater than the dimensions of the matrix.");
    }
    if (m1.numColumns <= col)
    {
        throw std::invalid_argument("The column being searched for is greater than the dimensions of the matrix.");
    }
    for (size_t i = m1.row_ptr.at(row); i < m1.row_ptr.at(row + 1); i++)
    {
        if (m1.col_ind.at(i) == col)
        {
            return m1.val.at(i);
        }
    }
    return 0;
}
/// @brief Converts a dense matrix to a compressed sparse row(CSR) matrix
/// @tparam T The type of the matrix
/// @param array The dense matrix to convert
/// @return The CSR matrix
template <typename T>
CSRMatrix<T> from_vector_CSR(vector<vector<T>> &array)
{
    CSRMatrix<T> returnMatrix;
    returnMatrix.numRows = array.size();
    returnMatrix.numColumns = array.at(0).size();
    returnMatrix.row_ptr.push_back(0);
    for (size_t i = 0; i < array.size(); i++)
    {
        for (size_t j = 0; j < array.at(i).size(); j++)
        {
            if (array.at(i).at(j) != 0)
            {
                returnMatrix.val.push_back(array.at(i).at(j));
                returnMatrix.col_ind.push_back(j);
            }
        }
        returnMatrix.row_ptr.push_back(returnMatrix.val.size());
    }
    return returnMatrix;
}

/// @brief Prints out a compressed sparse row(CSR) matrix to cout
/// @tparam T The type of the matrix
/// @param m1 The matrix too print out
template <typename T>
void print_matrixCSR(CSRMatrix<T> m1)
{
    for (size_t i = 0; i < m1.numRows; i++)
    {
        for (size_t j = 0; j < m1.numColumns; j++)
        {
            cout << get_matrixCSR(m1, i, j) << " ";
        }
        cout << endl;
    }
}
/// @brief Adds two compressed spares row(CSR) matrices together
/// @exception The two matrixes must have the same dimensions
/// @tparam T The type of both matrixes
/// @param m1 The first matrix too add
/// @param m2 The second matrix too add
/// @return m1+m2
template <typename T>
CSRMatrix<T> add_matrixCSR(CSRMatrix<T> m1, CSRMatrix<T> m2)
{
    if (m1.numRows != m2.numRows)
    {
        throw std::invalid_argument("The number of rows in the first matrix must match the number of rows in the second matrix.");
    }
    if (m1.numColumns != m2.numColumns)
    {
        throw std::invalid_argument("The number of columns in the first matrix must match the number of columns in the second matrix.");
    }
    CSRMatrix<T> returnMatrix;
    returnMatrix.numRows = m1.numRows;
    returnMatrix.numColumns = m1.numColumns;
    returnMatrix.row_ptr.push_back(0);

    for (size_t i = 0; i < m1.numRows; i++)
    {
        size_t a1 = m1.row_ptr.at(i);
        size_t b1 = m1.row_ptr.at(i + 1);
        size_t a2 = m2.row_ptr.at(i);
        size_t b2 = m2.row_ptr.at(i + 1);
        while (a1 < b1 && a2 < b2)
        {
            if (m1.col_ind.at(a1) < m2.col_ind.at(a2))
            {
                returnMatrix.val.push_back(m1.val.at(a1));
                returnMatrix.col_ind.push_back(m1.col_ind.at(a1));
                a1++;
            }
            else if (m1.col_ind.at(a1) > m2.col_ind.at(a2))
            {
                returnMatrix.val.push_back(m2.val.at(a2));
                returnMatrix.col_ind.push_back(m2.col_ind.at(a2));
                a2++;
            }
            else if (m1.col_ind.at(a1) == m2.col_ind.at(a2))
            {
                T value = m1.val.at(a1) + m2.val.at(a2);
                if (value != 0)
                {
                    returnMatrix.val.push_back(value);
                    returnMatrix.col_ind.push_back(m1.col_ind.at(a1));
                }
                a1++;
                a2++;
            }
        }
        while (a1 < b1)
        {
            returnMatrix.val.push_back(m1.val.at(a1));
            returnMatrix.col_ind.push_back(m1.col_ind.at(a1));
            a1++;
        }
        while (a2 < b2)
        {
            returnMatrix.val.push_back(m2.val.at(a2));
            returnMatrix.col_ind.push_back(m2.col_ind.at(a2));
            a2++;
        }
        returnMatrix.row_ptr.push_back(returnMatrix.val.size());
    }
    return returnMatrix;
}

/// @brief Transposes a compressed sparse row(CSR) matrix
/// @tparam T The type of the matrix
/// @param m1 The CSR matrix to transpose
/// @return m1 transposed
/// We first check make a vector row_count with 0s. This will keep track of how many elements
/// are at each column (or each row for the transposed matrix). We compute this by iterating through
/// the elements and adding the count at the respective row. Then, each element row_ptr[i] of the new
/// matrix row_ptr is simply the sum of the row_count up until i. Later on, we clear up row_count
/// back to 0s and again iterate through the elements and update row_count accordingly. In this case,
/// for an element in row i and overall index j in the original matrix we can compute its new overall
/// index in the transposed matrix as follows: Take the number of elements in the columns behind it:
/// returnMatrix.row_ptr[m1.col_ind[j]]. Sum it with the elements in this column so far:
/// row_count[m1.col_ind[j]]. And obtain the index. Then add this element to such index in the new
/// matrix: returnMatrix.val[index]=m1.val[j] as well as its column index: returnMatrix.col_ind[index]=i.
template <typename T>
CSRMatrix<T> transpose_matrixCSR(CSRMatrix<T> m1)
{
    CSRMatrix<T> returnMatrix;
    returnMatrix.numRows = m1.numColumns;
    returnMatrix.numColumns = m1.numRows;
    returnMatrix.row_ptr.push_back(0);
    vector<size_t> row_count(m1.numColumns, 0);
    returnMatrix.val = vector<T>(m1.val.size());
    returnMatrix.col_ind = vector<size_t>(m1.col_ind.size());
    for (size_t i = 0; i < m1.numRows; i++)
    {
        for (size_t j = m1.row_ptr.at(i); j < m1.row_ptr.at(i + 1); j++)
        {
            row_count.at(m1.col_ind.at(j))++;
        }
    }
    for (size_t i = 0; i < m1.numColumns; i++)
    {
        returnMatrix.row_ptr.push_back(returnMatrix.row_ptr.at(i) + row_count.at(i));
    }
    row_count = vector<size_t>(m1.numColumns, 0);
    for (size_t i = 0; i < m1.numRows; i++)
    {
        for (size_t j = m1.row_ptr.at(i); j < m1.row_ptr.at(i + 1); j++)
        {
            size_t index = returnMatrix.row_ptr.at(m1.col_ind.at(j)) + row_count.at(m1.col_ind.at(j));
            returnMatrix.val.at(index) = m1.val.at(j);
            returnMatrix.col_ind.at(index) = i;
            row_count.at(m1.col_ind.at(j))++;
        }
    }
    return returnMatrix;
}

/// @brief Multiplies two compressed sparse row(CSR) matrixes
/// @exception The number of columns in m1 must equal the number of rows in m2
/// @tparam T The type of the matrixes
/// @param m1 The first CSR matrix to multiply
/// @param m2 The second CSR matrix to multiply
/// @return The dot product of m1 and m2
template <typename T>
CSRMatrix<T> multiply_matrixCSR(CSRMatrix<T> m1, CSRMatrix<T> m2)
{
    if (m1.numColumns != m2.numRows)
    {
        throw std::invalid_argument("The number of columns in the first matrix must match the number of rows in the second matrix.");
    }
    CSRMatrix<T> returnMatrix;
    returnMatrix.numRows = m1.numRows;
    returnMatrix.numColumns = m2.numColumns;
    returnMatrix.row_ptr.push_back(0);
    CSRMatrix<T> m2t = transpose_matrixCSR(m2);
    for (size_t i = 0; i < m1.numRows; i++)
    {
        for (size_t j = 0; j < m2t.numRows; j++)
        {
            T sum = 0;
            size_t a1 = m1.row_ptr.at(i);
            size_t b1 = m1.row_ptr.at(i + 1);
            size_t a2 = m2t.row_ptr.at(j);
            size_t b2 = m2t.row_ptr.at(j + 1);
            while (a1 < b1 && a2 < b2)
            {
                if (m1.col_ind.at(a1) < m2t.col_ind.at(a2))
                {
                    a1++;
                }
                else if (m1.col_ind.at(a1) > m2t.col_ind.at(a2))
                {
                    a2++;
                }
                else if (m1.col_ind.at(a1) == m2t.col_ind.at(a2))
                {
                    sum += m1.val.at(a1) * m2t.val.at(a2);
                    a1++;
                    a2++;
                }
            }
            if (sum != 0)
            {
                returnMatrix.val.push_back(sum);
                returnMatrix.col_ind.push_back(j);
            }
        }
        returnMatrix.row_ptr.push_back(returnMatrix.val.size());
    }
    return returnMatrix;
}

/// @brief Subtract two compressed sparse row(CSR) matrixes
/// @exception The number of columns in m1 must equal the number of rows in m2
/// @tparam T The type of the matrixes
/// @param m1 The first CSR matrix to subtract
/// @param m2 The second CSR matrix to subtract
/// @return The difference of m1 and m2
template <typename T>
CSRMatrix<T> subtract_matrixCSR(CSRMatrix<T> m1, CSRMatrix<T> m2)
{
    if (m1.numRows != m2.numRows)
    {
        throw std::invalid_argument("The number of rows in the first matrix must match the number of rows in the second matrix.");
    }
    if (m1.numColumns != m2.numColumns)
    {
        throw std::invalid_argument("The number of columns in the first matrix must match the number of columns in the second matrix.");
    }
    CSRMatrix<T> returnMatrix;
    returnMatrix.numRows = m1.numRows;
    returnMatrix.numColumns = m1.numColumns;
    returnMatrix.row_ptr.push_back(0);

    for (size_t i = 0; i < m1.numRows; i++)
    {
        size_t a1 = m1.row_ptr.at(i);
        size_t b1 = m1.row_ptr.at(i + 1);
        size_t a2 = m2.row_ptr.at(i);
        size_t b2 = m2.row_ptr.at(i + 1);
        while (a1 < b1 && a2 < b2)
        {
            if (m1.col_ind.at(a1) < m2.col_ind.at(a2))
            {
                returnMatrix.val.push_back(m1.val.at(a1));
                returnMatrix.col_ind.push_back(m1.col_ind.at(a1));
                a1++;
            }
            else if (m1.col_ind.at(a1) > m2.col_ind.at(a2))
            {
                returnMatrix.val.push_back(-m2.val.at(a2));
                returnMatrix.col_ind.push_back(m2.col_ind.at(a2));
                a2++;
            }
            else if (m1.col_ind.at(a1) == m2.col_ind.at(a2))
            {
                T value = m1.val.at(a1) - m2.val.at(a2);
                if (value != 0)
                {
                    returnMatrix.val.push_back(value);
                    returnMatrix.col_ind.push_back(m1.col_ind.at(a1));
                }
                a1++;
                a2++;
            }
        }
        while (a1 < b1)
        {
            returnMatrix.val.push_back(m1.val.at(a1));
            returnMatrix.col_ind.push_back(m1.col_ind.at(a1));
            a1++;
        }
        while (a2 < b2)
        {
            returnMatrix.val.push_back(-m2.val.at(a2));
            returnMatrix.col_ind.push_back(m2.col_ind.at(a2));
            a2++;
        }
        returnMatrix.row_ptr.push_back(returnMatrix.col_ind.size());
    }
    return returnMatrix;
}

/// @brief Scalar multiply a compressed sparse row(CSR) matrix
/// @tparam T The type of the matrix
/// @param m The CSR matrix to scalar multiply
/// @param scalar The scalar to multiply the matrix by
/// @return The scalar multiplied matrix
template <typename T>
CSRMatrix<T> scalar_multiply_CSR(CSRMatrix<T> m, T scalar)
{

    CSRMatrix<T> result;
    result.numRows = m.numRows;
    result.numColumns = m.numColumns;
    result.val = vector<T>(m.val.size());
    result.col_ind = vector<size_t>(m.col_ind.begin(), m.col_ind.end());
    result.row_ptr = vector<size_t>(m.row_ptr.begin(), m.row_ptr.end());
    for (size_t i = 0; i < m.val.size(); i++)
    {
        result.val.at(i) = m.val.at(i) * scalar;
    }
    return result;
}

/// @brief Find the min value in a compressed sparse row(CSR) matrix
/// @tparam T The type of the matrix
/// @param m The CSR matrix to find the min value of
/// @return The min value in the matrix

template <typename T>
T find_min_CSR(CSRMatrix<T> matrix)
{

    T min_value = matrix.val[0];
    for (T val : matrix.val)
    {
        if (val < min_value)
        {
            min_value = val;
        }
    }
    return min_value;
}

/// @brief Find the max value in a compressed sparse row(CSR) matrix
/// @tparam T The type of the matrix
/// @param m The CSR matrix to find the max value of
/// @return The max value in the matrix
template <typename T>
T find_max_CSR(CSRMatrix<T> matrix)
{
    T max_value = matrix.val[0];
    for (T val : matrix.val)
    {
        if (val > max_value)
        {
            max_value = val;
        }
    }
    return max_value;
}

// IN CORRECT BECAUSE .mtx is not sorted by rows
// /// @brief Creates CSR matrix for .mtx file
// /// @tparam T the type of matrix
// /// @param fileName the name of the file to import
// /// @return a new CSR matrix from the filename
// template <typename T>
// CSRMatrix<T> load_fileCSR(string fileName)
// {
//     std::ifstream file(fileName);
//     int num_row = 0, num_col = 0, num_lines = 0;

//     // Ignore comments headers
//     while (file.peek() == '%')
//         file.ignore(2048, '\n');

//     // Read number of rows, columns, and non-zero values
//     file >> num_row >> num_col >> num_lines;

//     CSRMatrix<T> returnMatrix;
//     returnMatrix.numRows = num_row;
//     returnMatrix.numColumns = num_col;
//     returnMatrix.row_ptr.push_back(0);
//     T data;
//     int row, col;
//     //Need to sort rows to ensure in row ascending order
//     file >> row >> col >> data;
//     row = row - 1;
//     for (int i = 0; i < num_row; i++)
//     {
//         while (num_lines > 0 && row == i)
//         {
//             returnMatrix.val.push_back(data);
//             returnMatrix.col_ind.push_back(col--);
//             file >> row >> col >> data;
//             row--;
//             num_lines--;
//         }
//         // this happens every time
//         returnMatrix.row_ptr.push_back(returnMatrix.val.size());
//     }

//     file.close();

//     return returnMatrix;
// }

//reincluding dense load, importing functions.cc causes redefination errors
/// @brief Creates dense matrix for .mtx file
/// @tparam T the type of matrix
/// @param fileName the name of the file to import
/// @return a new dense matrix from the filename
template <typename T>
std::vector<std::vector<T>> load_fileMatrix_internal(const std::string& fileName) {
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
//SLOW FILE READ BUT FINE FOR SMALL MATRICES
/// @brief Creates CSR matrix for .mtx file
/// @tparam T the type of matrix
/// @param fileName the name of the file to import
/// @return a new CSR matrix from the filename
template <typename T>
CSRMatrix<T> load_fileCSR(string fileName)
{
    vector<vector<T>> dense_result = load_fileMatrix_internal<T>(fileName);
    return from_vector_CSR<T>(dense_result);
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
template <typename T>
bool diagonally_dominant(CSRMatrix<T> m1) {
    for (size_t i = 0; i < m1.numRows; ++i) {
        size_t a1 = m1.row_ptr.at(i);
        size_t b1 = m1.row_ptr.at(i + 1);
        T sum = 0.0;
        T diagonal = 0.0;
        while(a1 < b1){
            if(m1.col_ind[a1] == i){
                diagonal = std::abs(m1.val[a1]);
            }else{
                sum += std::abs(m1.val[a1]);
            }
            a1++;
        }
        if (diagonal < sum) {
            return false;
        }
    }
    
    return true;
}
#include <cmath>

/**
 * @brief The Jacobi Method is an iterative method for determining the solutions of a strictly
 * diagonally dominant matrix A. Through each iteration, the values of x[i] are approximated through
 * the formula x[i] = B[i]
 * 
 * @param denseMatrix 
 * @param B 
 * @param tol - the tolerance for convergence
 * @param iterations - the maximum number of iterations to perform
 */
template <typename T>
std::vector<T> jacobi_method_CSR(CSRMatrix<T> m1, std::vector<T> B, const double tol,int maxIterations) {
    // if (diagonally_dominant(m1) == false) {
    //     throw std::invalid_argument("Input matrix is not diagonally dominant");
    // }
    std::vector<T> xValues(B.size(), 0.0);
    std::vector<T> approxValues(B.size(), 0.0);
    int iterations = 0;
    double diff = tol + 1.0;
    while (iterations < maxIterations && diff > tol) {
        for (size_t i = 0; i < m1.numRows; ++i) {
            size_t a1 = m1.row_ptr.at(i);
            size_t b1 = m1.row_ptr.at(i + 1);
            T sum = 0.0;
            T diagonal = 0.0;
            while(a1 < b1){
                if(m1.col_ind[a1] == i){
                    diagonal = m1.val[a1];
                }else{
                    sum += m1.val[a1] * xValues[m1.col_ind[a1]];
                }
                a1++;
            }
            //no devide by zero error becuase of diagonally dominant check
            approxValues[i] = (B[i] - sum) / diagonal;
            
        }
         // Calculate the norm of the difference between x and x_new
        diff = 0.0;
        for (size_t i = 0; i < m1.numRows; i++)
        {
            double abs_diff = std::abs(approxValues[i] - xValues[i]);
            if (abs_diff > diff)
            {
                diff = abs_diff;
            }
        }
        xValues = approxValues;
        iterations++;
}
    cerr << "Jacobi Sparse Itertions: "<< iterations <<endl;

return approxValues;
}

/**
 * @brief CSR Gauss-Seidel Method. Similar to the Jacobi method however we update the X vector directly
 * instead
 * 
 * @tparam T 
 * @param m1 
 * @param B 
 * @param tol 
 * @param maxIterations 
 * @return std::vector<T> 
 */
template <typename T>
std::vector<T> gauss_sidel_CSR(CSRMatrix<T> m1, std::vector<T> B, const double tol,int maxIterations) {
    // if (diagonally_dominant(m1) == false) {
    //     throw std::invalid_argument("Input matrix is not diagonally dominant");
    // }
    std::vector<T> xValues(B.size(), 0.0);
    std::vector<T> approxValues(B.size(), 0.0);
    int iterations = 0;
    double diff = tol + 1.0;
    while (iterations < maxIterations && diff > tol) {
        for (size_t i = 0; i < m1.numRows; ++i) {
            size_t a1 = m1.row_ptr.at(i);
            size_t b1 = m1.row_ptr.at(i + 1);
            T sum = 0.0;
            T diagonal = 0.0;
            while(a1 < b1){
                if(m1.col_ind[a1] == i){
                    diagonal = m1.val[a1];
                }else{
                    sum += m1.val[a1] * xValues[m1.col_ind[a1]];
                }
                a1++;
            }
            //no divide by zero error becuase of diagonally dominant check
            approxValues[i] = (B[i] - sum) / diagonal;
            
        }
        diff = 0.0;
        for (int i = 0; i < m1.numRows; i++)
        {
            T abs_diff = std::abs(approxValues[i] - xValues[i]);
            if (abs_diff > diff)
            {
                diff = abs_diff;
            }
        }
        xValues = approxValues;
        iterations++;
    }
    cerr << "Gauss Sidel Sparse Itertions: "<< iterations <<endl;
    return xValues;
}

template <typename T>
std::vector<T> ssor_iteration_CSR(CSRMatrix<T> A,
                                  const std::vector<T> &b,
                                  const T tol,
                                  const int max_iter,
                                  const T omega)
{
    const int n = A.numRows;
    std::vector<T> x(n, 0.0);
    std::vector<T> x_new(n, 0.0);

    int iter = 0;
    T diff = tol + 1.0;

    while (iter < max_iter && diff > tol)
    {
        // Normally one needs to solve M(x_new - x) = b - Ax using lu_solve where M = (D + omega*L)^-1 * (D + omega*U),
        // L is the strict lower triangular part of A, U is the strict upper triangular part of A, and D is the diagonal of A.
        // However, according to https://en.wikipedia.org/wiki/Successive_over-relaxation, we can use forward substitution:
        for (int i = 0; i < n; i++)
        {
            T sum = 0.0, diag;
            for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; k++) {
                int j = A.col_ind[k];
                if (j < i) {
                    sum += A.val[k] * x_new[j];
                } else if (j > i) {
                    sum += A.val[k] * x[j];
                } else {
                    diag = A.val[k];
                }
            }
            x_new[i] = (1.0 - omega) * x[i] + (omega / diag) * (b[i] - sum);
        }
        // Compute the difference between the new and old iterates
        diff = 0.0;
        for (int i = 0; i < n; i++)
        {
            T abs_diff = std::abs(x_new[i] - x[i]);
            if (abs_diff > diff)
            {
                diff = abs_diff;
            }
        }
        x = x_new;
        iter++;
    }
    cerr << "SSOR Sparse Itertions: "<< iter <<endl;
    return x;
}


/**
 * @brief 
 * 
 * @tparam T 
 * @param m1 
 */
template <typename T> 
    void lu_decomposition_CSR(CSRMatrix<T> m1) {
        int n = m1.row_ptr.size() - 1;
        CSRMatrix<T> L;
        L.val = m1.vals;
        L.row_ptr = m1.row_ptr;
        L.col_ind = m1.col_ind;
        CSRMatrix<T> U;
        U.val = m1.val;
        U.row_ptr = m1.row_ptr;
        U.col_ind = m1.col_ind;

        for (int i = 0; i < n; ++i) {
            L.val.at(L.row_ptr[i]) = 1.0;

            
        }


    }

// int main() {
//     vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};
//     vector<vector<int> > array2 = {{0,2,3},{0,-5,0},{7,8,0}};
//     CSRMatrix<int> matrix = from_vector(array);
//     CSRMatrix<int> matrix2 = from_vector(array2);

//     CSRMatrix<int> matrix3 = add_matrixCSR(matrix, matrix2);
//     CSRMatrix<int> matrix4 = scalar_multiply_CSR(matrix, 2);
//     CSRMatrix<int> matrix5 = scalar_multiply_CSR(matrix2, 2);
//     CSRMatrix<int> matrix6 = subtract_matrixCSR(matrix, matrix2);

//     cout << "Matrix 1" << endl;
//     print_matrixCSR(matrix);
//     cout << "Matrix 2" << endl;
//     print_matrixCSR(matrix2);
//     cout << "Matrix 3" << endl;
//     print_matrixCSR(matrix3);
//     cout << "Matrix 4" << endl;
//     print_matrixCSR(matrix4);
//     cout << "Matrix 5" << endl;
//     print_matrixCSR(matrix5);
//     cout << "Matrix 6" << endl;
//     print_matrixCSR(matrix6);
//     cout << find_min_CSR(matrix) << endl;
//     cout << find_max_CSR(matrix) << endl;

//     return 0;
// }
#include <iostream>
#include <vector>
#include <cassert>

#include <chrono>
// class timer {
// public:
//     std::chrono::time_point<std::chrono::high_resolution_clock> lastTime;
//     timer() : lastTime(std::chrono::high_resolution_clock::now()) {}
//     inline double elapsed() {
//         std::chrono::time_point<std::chrono::high_resolution_clock> thisTime=std::chrono::high_resolution_clock::now();
//         double deltaTime = std::chrono::duration<double>(thisTime-lastTime).count();
//         lastTime = thisTime;
//         return deltaTime;
//     }
// };
// #include <iomanip>

// void printVector(vector<double> &B) {
//     // Set the width of each output element to 8 characters
//     const int width = 2;
//         // Loop over each element in the row
//         for (size_t j = 0; j < B.size(); j++) {
//             // Set the output width for the element
//             cout << setw(width) << B[j];
//         }
//         // End the row with a newline
//         cout << endl;
// }

void compareVector(vector<double> &A, vector<double> &B, double error) {
    if (A.size () != B.size()){
        cerr<< "NOT SAME SIZE " <<endl;
    }
        // Loop over each element in the row
        for (size_t i = 0; i < B.size(); i++) {
            if(A[i] - B[i] > error || B[i] - A[i] > error){
                cerr<< "NOT THE SAME " <<endl;
            }
        }
}

int main() {
    std::string to_load = "bcsstk10.mtx";
    CSRMatrix<double> A_CSR = load_fileCSR<double>(to_load);
    vector<vector<double>> A = load_fileMatrix<double>(to_load);
    vector<double> B(A_CSR.numRows, 0);
    B[0] = 1.0;
    timer stopwatch;

    // stopwatch.elapsed();
    // gcr(A,B,1e-6,1000);
    // cerr << "GCR TIME RESULTS SECONDS: "<< stopwatch.elapsed() <<endl <<endl;

    stopwatch.elapsed();
    jacobi_iteration(A,B,1e-16,1000);
    cerr << "Jacobi TIME RESULTS SECONDS: "<< stopwatch.elapsed() <<endl;
    stopwatch.elapsed();
    jacobi_method_CSR<double>(A_CSR,B,1e-16,1000);
    cerr << "Jacobi CSR TIME RESULTS SECONDS: "<< stopwatch.elapsed() <<endl <<endl;


    stopwatch.elapsed();
    gauss_seidel(A,B,1e-16,1000);
    cerr << "Gauss_seidel TIME RESULTS SECONDS: "<< stopwatch.elapsed() <<endl;
    stopwatch.elapsed();
    gauss_sidel_CSR<double>(A_CSR,B,1e-16,1000);
    cerr << "Gauss Seidel CSR TIME RESULTS SECONDS: "<< stopwatch.elapsed() <<endl <<endl;

    stopwatch.elapsed();
    ssor_iteration(A,B,1e-16,1000,0.5);
    cerr << "SSOR w = 0.5 TIME RESULTS SECONDS: "<< stopwatch.elapsed() <<endl;
    stopwatch.elapsed();
    ssor_iteration_CSR<double>(A_CSR,B,1e-16,1000,0.5);
    cerr << "SSOR w = 0.5 CSR TIME RESULTS SECONDS: "<< stopwatch.elapsed() <<endl <<endl;

    // stopwatch.elapsed();
    // ssor_iteration(A,B,1e-16,1000,1.5);
    // cerr << "SSOR w = 1.5 TIME RESULTS SECONDS: "<< stopwatch.elapsed() <<endl;
    // stopwatch.elapsed();
    // ssor_iteration_CSR<double>(A_CSR,B,1e-16,1000,1.5);
    // cerr << "SSOR w = 1.5 CSR TIME RESULTS SECONDS: "<< stopwatch.elapsed() <<endl <<endl;

    // stopwatch.elapsed();
    // ssor_iteration(A,B,1e-16,1000,2.0);
    // cerr << "SSOR w = 2.0 TIME RESULTS SECONDS: "<< stopwatch.elapsed() <<endl;
    // stopwatch.elapsed();
    // ssor_iteration_CSR<double>(A_CSR,B,1e-16,1000,2.0);
    // cerr << "SSOR w = 2.0 CSR TIME RESULTS SECONDS: "<< stopwatch.elapsed() <<endl <<endl;

    stopwatch.elapsed();
    gaussian_elimination(A,B);
    cerr << "Gauss Elimination TIME RESULTS SECONDS: "<< stopwatch.elapsed() <<endl;
    return 0;
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
//     //std::vector<int> size = {5000,10000,15000,20000,25000,30000,35000,40000};
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
//         // CSRMatrix<double> A_CSR = from_vector_CSR<double>(A);
//         // stopwatch.elapsed();
//         // jacobi_method_CSR<double>(A_CSR,B,1e-6,1000);
//         // jacobi.push_back(stopwatch.elapsed());
//         // stopwatch.elapsed();
//         // ssor_iteration_CSR<double>(A_CSR,B,1e-6,1000,0.5);
//         // omega_5.push_back(stopwatch.elapsed());
//         // stopwatch.elapsed();
//         // gauss_sidel_CSR<double>(A_CSR,B,1e-6,1000);
//         // gauss_seidel_time.push_back(stopwatch.elapsed());


//         CSRMatrix<double> A_CSR = from_vector_CSR<double>(A);
//         jacobi = jacobi_method_CSR<double>(A_CSR,B,1e-16,1000);
//     }
//     cerr<< "Jacobi CSR: ";
//     for(auto val : jacobi){
//         cerr<< val << ",";
//     }
//     cerr<< endl;
//     // cerr<< "SSOR CSR w = 0.5: ";
//     // for(auto val : omega_5){
//     //     cerr<< val << ",";
//     // }
//     // cerr<< endl;
//     // cerr<< "Gauss Seidel CSR: ";
//     // for(auto val : gauss_seidel_time){
//     //     cerr<< val << ",";
//     // }
//     // cerr<< endl;

//   return 0;
// }