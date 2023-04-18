#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

template <typename T>
class CSCMatrix
{
public:
    size_t numRows, numColumns;
    vector<T> val;
    vector<size_t> row_ind;
    vector<size_t> col_ptr;
};
// loadfile
// savefile
/// @brief Gets a value from the compressed sparse column(CSC) matrix
/// @exception The search row and col must be less than dimensions of m1
/// @tparam T The type of the matrix
/// @param m1 The CSC matrix to get the value from
/// @param row The row of the value to get
/// @param col The column of the value to get
/// @return That value stored at row,column
template <typename T>
T get_matrixCSC(CSCMatrix<T> m1, size_t row, size_t col)
{
    if (m1.numRows <= row)
    {
        // throw error
        cout << "Number of rows exceeded: " << row << endl;
    }
    if (m1.numColumns <= col)
    {
        // throw error
        cout << "Number of columns exceeded: " << col << endl;
    }
    // return -1;
    for (size_t i = m1.col_ptr.at(col); i < m1.col_ptr.at(col + 1); i++)
    {
        if (m1.row_ind.at(i) == row)
        {
            return m1.val.at(i);
        }
    }
    return 0;
}

/// @brief Converts a dense matrix to a compressed sparse column(CSC) matrix
/// @tparam T The type of the matrix
/// @param array  The dense matrix to convert
/// @return The CSC matrix
template <typename T>
CSCMatrix<T> from_vector(vector<vector<T>> &array)
{
    CSCMatrix<T> returnMatrix;
    returnMatrix.numRows = array.size();
    returnMatrix.numColumns = array.at(0).size();
    returnMatrix.col_ptr.push_back(0);
    for (size_t i = 0; i < array.at(0).size(); i++)
    {
        for (size_t j = 0; j < array.size(); j++)
        {
            if (array.at(j).at(i) != 0)
            {
                returnMatrix.val.push_back(array.at(j).at(i));
                returnMatrix.row_ind.push_back(j);
            }
        }
        returnMatrix.col_ptr.push_back(returnMatrix.val.size());
    }
    return returnMatrix;
}

/// @brief Converts a compressed sparse column(CSC) matrix to a dense matrix
/// @tparam T The type of the matrix
/// @param m1 The CSC matrix to convert
/// @return The dense matrix
template <typename T>
void print_matrixCSC(CSCMatrix<T> m1)
{
    for (size_t i = 0; i < m1.numRows; i++)
    {
        for (size_t j = 0; j < m1.numColumns; j++)
        {
            cout << get_matrixCSC(m1, i, j) << " ";
        }
        cout << endl;
    }
}

/// @brief Adds two compressed sparse column(CSC) matrices
/// @tparam T The type of the matrix
/// @param m1 The first CSC matrix to add
/// @param m2 The second CSC matrix to add
/// @return The sum of the two matrices
template <typename T>
CSCMatrix<T> add_matrixCSC(CSCMatrix<T> m1, CSCMatrix<T> m2)
{
    if (m1.numRows != m2.numRows)
    {
        // throw error
        cout << "Number of rows not equal: " << m1.numRows << " " << m2.numRows << endl;
    }
    if (m1.numColumns != m2.numColumns)
    {
        // throw error
        cout << "Number of columns not equal: " << m1.numColumns << " " << m2.numColumns << endl;
    }
    CSCMatrix<T> returnMatrix;
    returnMatrix.numRows = m1.numRows;
    returnMatrix.numColumns = m1.numColumns;
    returnMatrix.col_ptr.push_back(0);

    for (size_t i = 0; i < m1.numColumns; i++)
    {
        size_t a1 = m1.col_ptr.at(i);
        size_t b1 = m1.col_ptr.at(i + 1);
        size_t a2 = m2.col_ptr.at(i);
        size_t b2 = m2.col_ptr.at(i + 1);
        while (a1 < b1 && a2 < b2)
        {
            if (m1.row_ind.at(a1) < m2.row_ind.at(a2))
            {
                returnMatrix.val.push_back(m1.val.at(a1));
                returnMatrix.row_ind.push_back(m1.row_ind.at(a1));
                a1++;
            }
            else if (m1.row_ind.at(a1) > m2.row_ind.at(a2))
            {
                returnMatrix.val.push_back(m2.val.at(a2));
                returnMatrix.row_ind.push_back(m2.row_ind.at(a2));
                a2++;
            }
            else if (m1.row_ind.at(a1) == m2.row_ind.at(a2))
            {
                T value = m1.val.at(a1) + m2.val.at(a2);
                if (value != 0)
                {
                    returnMatrix.val.push_back(value);
                    returnMatrix.row_ind.push_back(m1.row_ind.at(a1));
                }
                a1++;
                a2++;
            }
        }
        while (a1 < b1)
        {
            returnMatrix.val.push_back(m1.val.at(a1));
            returnMatrix.row_ind.push_back(m1.row_ind.at(a1));
            a1++;
        }
        while (a2 < b2)
        {
            returnMatrix.val.push_back(m2.val.at(a2));
            returnMatrix.row_ind.push_back(m2.row_ind.at(a2));
            a2++;
        }
        returnMatrix.col_ptr.push_back(returnMatrix.val.size());
    }
    return returnMatrix;
}

/// @brief Transposes a compressed sparse column(CSC) matrix
/// @tparam T  The type of the matrix
/// @param m1  The CSC matrix to transpose
/// @return     The transposed CSC matrix
template <typename T>
CSCMatrix<T> transpose_matrixCSC(CSCMatrix<T> m1)
{
    CSCMatrix<T> returnMatrix;
    returnMatrix.numRows = m1.numColumns;
    returnMatrix.numColumns = m1.numRows;
    returnMatrix.col_ptr.push_back(0);
    vector<size_t> row_count(m1.numRows, 0);
    returnMatrix.val = vector<T>(m1.val.size());
    returnMatrix.row_ind = vector<size_t>(m1.row_ind.size());
    for (size_t i = 0; i < m1.numColumns; i++)
    {
        for (size_t j = m1.col_ptr.at(i); j < m1.col_ptr.at(i + 1); j++)
        {
            row_count.at(m1.row_ind.at(j))++;
        }
    }
    for (size_t i = 0; i < m1.numRows; i++)
    {
        returnMatrix.col_ptr.push_back(returnMatrix.col_ptr.at(i) + row_count.at(i));
    }
    row_count = vector<size_t>(m1.numRows, 0);
    for (size_t i = 0; i < m1.numColumns; i++)
    {
        for (size_t j = m1.col_ptr.at(i); j < m1.col_ptr.at(i + 1); j++)
        {
            size_t index = returnMatrix.col_ptr.at(m1.row_ind.at(j)) + row_count.at(m1.row_ind.at(j));
            returnMatrix.val.at(index) = m1.val.at(j);
            returnMatrix.row_ind.at(index) = i;
            row_count.at(m1.row_ind.at(j))++;
        }
    }
    return returnMatrix;
}

/// @brief Multiplies two compressed sparse column(CSC) matrices
/// @tparam T The type of the matrix
/// @param m1 The first CSC matrix to multiply
/// @param m2 The second CSC matrix to multiply
/// @return The product of the two matrices
template <typename T>
CSCMatrix<T> multiply_matrixCSC(CSCMatrix<T> m1, CSCMatrix<T> m2)
{
    if (m1.numColumns != m2.numRows)
    {
        // throw error
        cout << "Number of columns of first matrix not equal to number of rows of second matrix: " << m1.numColumns << " " << m2.numRows << endl;
    }
    CSCMatrix<T> returnMatrix;
    returnMatrix.numRows = m1.numRows;
    returnMatrix.numColumns = m2.numColumns;
    returnMatrix.col_ptr.push_back(0);
    CSCMatrix<T> m1t = transpose_matrixCSC(m1);
    for (size_t i = 0; i < m1t.numColumns; i++)
    {
        for (size_t j = 0; j < m2.numColumns; j++)
        {
            T sum = 0;
            size_t a1 = m1t.col_ptr.at(i);
            size_t b1 = m1t.col_ptr.at(i + 1);
            size_t a2 = m2.col_ptr.at(j);
            size_t b2 = m2.col_ptr.at(j + 1);
            while (a1 < b1 && a2 < b2)
            {
                if (m1t.row_ind.at(a1) < m2.row_ind.at(a2))
                {
                    a1++;
                }
                else if (m1t.row_ind.at(a1) > m2.row_ind.at(a2))
                {
                    a2++;
                }
                else if (m1t.row_ind.at(a1) == m2.row_ind.at(a2))
                {
                    sum += m1t.val.at(a1) * m2.val.at(a2);
                    a1++;
                    a2++;
                }
            }
            if (sum != 0)
            {
                returnMatrix.val.push_back(sum);
                returnMatrix.row_ind.push_back(j);
            }
        }
        returnMatrix.col_ptr.push_back(returnMatrix.val.size());
    }
    return returnMatrix;
}

/// @brief Subtract two compressed sparse column(CSC) matrices
/// @tparam T The type of the matrix
/// @param m1 The first CSC matrix to add
/// @param m2 The second CSC matrix to add
/// @return The difference of the two matrices
template <typename T>
CSCMatrix<T> subtract_matrixCSC(CSCMatrix<T> m1, CSCMatrix<T> m2)
{
    if (m1.numRows != m2.numRows)
    {
        // throw error
        cout << "Number of rows not equal: " << m1.numRows << " " << m2.numRows << endl;
    }
    if (m1.numColumns != m2.numColumns)
    {
        // throw error
        cout << "Number of columns not equal: " << m1.numColumns << " " << m2.numColumns << endl;
    }
    CSCMatrix<T> returnMatrix;
    returnMatrix.numRows = m1.numRows;
    returnMatrix.numColumns = m1.numColumns;
    returnMatrix.col_ptr.push_back(0);

    for (size_t i = 0; i < m1.numColumns; i++)
    {
        size_t a1 = m1.col_ptr.at(i);
        size_t b1 = m1.col_ptr.at(i + 1);
        size_t a2 = m2.col_ptr.at(i);
        size_t b2 = m2.col_ptr.at(i + 1);
        while (a1 < b1 && a2 < b2)
        {
            if (m1.row_ind.at(a1) < m2.row_ind.at(a2))
            {
                returnMatrix.val.push_back(m1.val.at(a1));
                returnMatrix.row_ind.push_back(m1.row_ind.at(a1));
                a1++;
            }
            else if (m1.row_ind.at(a1) > m2.row_ind.at(a2))
            {
                returnMatrix.val.push_back(-m2.val.at(a2));
                returnMatrix.row_ind.push_back(m2.row_ind.at(a2));
                a2++;
            }
            else if (m1.row_ind.at(a1) == m2.row_ind.at(a2))
            {
                T value = m1.val.at(a1) - m2.val.at(a2);
                if (value != 0)
                {
                    returnMatrix.val.push_back(value);
                    returnMatrix.row_ind.push_back(m1.row_ind.at(a1));
                }
                a1++;
                a2++;
            }
        }
        while (a1 < b1)
        {
            returnMatrix.val.push_back(m1.val.at(a1));
            returnMatrix.row_ind.push_back(m1.row_ind.at(a1));
            a1++;
        }
        while (a2 < b2)
        {
            returnMatrix.val.push_back(-m2.val.at(a2));
            returnMatrix.row_ind.push_back(m2.row_ind.at(a2));
            a2++;
        }
        returnMatrix.col_ptr.push_back(returnMatrix.val.size());
    }
    return returnMatrix;
}

/// @brief Scalar multiply a compressed sparse column(CSC) matrix
/// @tparam T The type of the matrix
/// @param m The CSC matrix to scalar multiply
/// @param scalar The scalar to multiply the matrix by
/// @return The scalar multiplied matrix
template <typename T>
CSCMatrix<T> scalar_multiply_CSC(CSCMatrix<T> m, T scalar)
{

    CSCMatrix<T> result;
    result.numRows = m.numRows;
    result.numColumns = m.numColumns;
    result.val = vector<T>(m.val.size());
    result.row_ind = vector<size_t>(m.row_ind.begin(), m.row_ind.end());
    result.col_ptr = vector<size_t>(m.col_ptr.begin(), m.col_ptr.end());
    for (size_t i = 0; i < m.val.size(); i++)
    {
        result.val.at(i) = m.val.at(i) * scalar;
    }
    return result;
}

/// @brief Find the min value in a compressed sparse row(CSC) matrix
/// @tparam T The type of the matrix
/// @param m The CSC matrix to find the min value of
/// @return The min value in the matrix
template <typename T>
T find_min_CSC(CSCMatrix<T> matrix)
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

/// @brief Find the max value in a compressed sparse row(CSC) matrix
/// @tparam T The type of the matrix
/// @param m The CSC matrix to find the max value of
/// @return The max value in the matrix
template <typename T>
T find_max_CSC(CSCMatrix<T> matrix)
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

int main()
{
    vector<vector<int>> array = vector<vector<int>>(4, vector<int>(3));
    array[0][0] = 1;
    array[0][1] = 2;
    array[0][2] = 3;
    array[1][0] = 4;
    array[1][1] = 5;
    array[1][2] = 6;
    array[2][0] = 7;
    array[2][1] = 8;
    array[2][2] = 9;
    array[3][0] = 10;
    array[3][1] = 11;
    array[3][2] = 12;

    vector<vector<int>> array2 = vector<vector<int>>(4, vector<int>(3));
    array2[0][0] = 9;
    array2[0][1] = -4;
    array2[0][2] = 7;
    array2[1][0] = -6;
    array2[1][1] = -5;
    array2[1][2] = 4;
    array2[2][0] = 3;
    array2[2][1] = 2;
    array2[2][2] = 1;
    array2[3][0] = 10;
    array2[3][1] = 11;
    array2[3][2] = 12;

    CSCMatrix<int> m1 = from_vector(array);
    CSCMatrix<int> m2 = from_vector(array2);
    CSCMatrix<int> m3 = add_matrixCSC(m1, m2);
    print_matrixCSC(m3);
    m3 = transpose_matrixCSC(m3);
    print_matrixCSC(m3);
    // m3 = transpose_matrixCSC(m3);
    m3 = multiply_matrixCSC(m2, transpose_matrixCSC(m1));
    print_matrixCSC(m3);
    cout << find_min_CSC(m3) << endl;
    cout << find_max_CSC(m3) << endl;
    print_matrixCSC(scalar_multiply_CSC(m3, 2));

    return 0;
}