#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

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
CSRMatrix<T> from_vector(vector<vector<T>> &array)
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

/// @brief Creates CSR matrix for .mtx file
/// @tparam T the type of matrix
/// @param fileName the name of the file to import
/// @return a new CSR matrix from the filename
template<typename T>
    CSRMatrix<T> load_fileCSR(string fileName){
    std::ifstream file(fileName);
    int num_row = 0, num_col = 0, num_lines = 0;

    // Ignore comments headers
    while (file.peek() == '%') file.ignore(2048, '\n');

    // Read number of rows, columns, and non-zero values
    file >> num_row >> num_col >> num_lines;

    CSRMatrix<T> returnMatrix;
    returnMatrix.numRows = num_row;
    returnMatrix.numColumns = num_col;
    returnMatrix.row_ptr.push_back(0);
    T data;
    int row, col;
    file >> row >> col >> data;
    row = row-1;
    for(int i = 0; i <num_row;i++){
        while(num_lines > 0 && row == i){
            returnMatrix.val.push_back(data);
            returnMatrix.col_ind.push_back(col);
            file >> row >> col >> data;
            row--;
            num_lines--;
        }
        //this happens every time
        returnMatrix.row_ptr.push_back(returnMatrix.val.size());
    }

    file.close();

    return returnMatrix;
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