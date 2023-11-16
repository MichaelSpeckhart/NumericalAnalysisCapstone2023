#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
// #include "../functionsCSRParallel.cc"
#include "../functions.cc"
#include "../functionsCSR.cc"
#include "../functionsCSC.cc"
#include "fstream"
const int numWidth = 10;
const char separator = ' ';
// Basic Unit tests for CSR add, multiply, and transpose
// Use -d to time the tests

/*
Takes a CSRMatrix and a dense matrix and makes sure they have the same elements
*/
/// @brief Ensure compressed sparse row matrix is the same as expected dense matrix
/// @param mResult The CSR matrix to compare to
/// @param mCheck The dense matrix to compare to
void CHECKCSR(CSRMatrix<int> mResult, vector<vector<int>> mCheck)
{
    CHECK(mResult.numRows == mCheck.size());
    for (size_t i = 0; i < mResult.numRows; i++)
    {
        CHECK(mResult.numColumns == mCheck[i].size());
        for (size_t j = 0; j < mResult.numColumns; j++)
        {
            CHECK_MESSAGE(get_matrixCSR(mResult, i, j) == mCheck[i][j], "i = " << i << ", j = " << j);
        }
    }
}

void CHECKCSR(CSRMatrix<double> mResult, vector<vector<double>> mCheck)
{
    CHECK(mResult.numRows == mCheck.size());
    for (size_t i = 0; i < mResult.numRows; i++)
    {
        CHECK(mResult.numColumns == mCheck[i].size());
        for (size_t j = 0; j < mResult.numColumns; j++)
        {
            CHECK_MESSAGE(get_matrixCSR(mResult, i, j) == mCheck[i][j], "i = " << i << ", j = " << j);
        }
    }
}

void CHECKMATRIX(vector<vector<double>> mResult, vector<vector<double>> mCheck)
{
    CHECK(mResult.size() == mCheck.size());
    for (size_t i = 0; i < mResult.size(); i++)
    {
        CHECK(mResult[i].size() == mCheck[i].size());
        for (size_t j = 0; j < mResult[i].size(); j++)
        {
            CHECK_MESSAGE(mResult[i][j] == mCheck[i][j], "i = " << i << ", j = " << j);
        }
    }
}

void CHECK_MATRIX_EQ(vector<vector<double>> &mResult, vector<vector<double>> &mCheck, double tol)
{
    CHECK(mResult.size() == mCheck.size());
    for (size_t i = 0; i < mResult.size(); i++)
    {
        CHECK(mResult[i].size() == mCheck[i].size());
        for (size_t j = 0; j < mResult[i].size(); j++)
        {
            CHECK_MESSAGE(abs(mResult[i][j] - mCheck[i][j]) < tol, "i = " << i << ", j = " << j);
        }
    }
}

void CHECK_VECTOR_EQ(vector<double> &mResult, vector<double> &mCheck, double tol)
{
    CHECK(mResult.size() == mCheck.size());
    for (size_t i = 0; i < mResult.size(); i++)
    {
        CHECK_MESSAGE(abs(mResult[i] - mCheck[i]) < tol, "i = " << i);
    }
}

template<typename T> void printElement(T t, const int& width)
{
    cout << left << setw(width) << setfill(separator) << t;
}

void PRINT_MATRIX(vector<vector<double>> &mResult)
{
    for (size_t i = 0; i < mResult.size(); i++) {
        for (size_t j = 0; j < mResult[0].size(); j++) {
            printElement(mResult[i][j], numWidth);  
        }
        cout << endl;
    }
}
TEST_CASE("testing CSR Add")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};

    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}};

    CSRMatrix<int> m1 = from_vector_CSR<int>(array);
    // check m1
    CHECKCSR(m1, array);

    CSRMatrix<int> m2 = from_vector_CSR<int>(array2);
    // check m2
    CHECKCSR(m2, array2);

    CSRMatrix<int> m3 = add_matrixCSR<int>(m1, m2);
    vector<vector<int>> addResultExpected = {{1, 2, 3}, {4, 0, 6}, {7, 16, 9}};
    // check the addition
    CHECKCSR(m3, addResultExpected);
}

TEST_CASE("CSR Add Exceptions")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}};

    CSRMatrix<int> m1 = from_vector_CSR<int>(array);
    CSRMatrix<int> m2 = from_vector_CSR<int>(array2);

    CHECK_THROWS_WITH_AS(add_matrixCSR<int>(m1, m2), "The number of rows in the first matrix must match the number of rows in the second matrix.", std::exception);

    vector<vector<int>> array3 = {{1, 0}, {4, 5}, {0, 8}};
    CSRMatrix<int> m3 = from_vector_CSR<int>(array3);

    CHECK_THROWS_WITH_AS(add_matrixCSR<int>(m1, m3), "The number of columns in the first matrix must match the number of columns in the second matrix.", std::exception);
}

TEST_CASE("testing CSR transpose")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};

    CSRMatrix<int> m1 = from_vector_CSR<int>(array);

    CSRMatrix<int> m2 = transpose_matrixCSR<int>(m1);
    vector<vector<int>> transposeResultExpected = {{1, 4, 0}, {0, 5, 8}, {0, 6, 9}};
    // check the transpose
    CHECKCSR(m2, transposeResultExpected);
}

TEST_CASE("testing CSR multiply")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}};

    CSRMatrix<int> m1 = from_vector_CSR<int>(array);
    CSRMatrix<int> m2 = from_vector_CSR<int>(array2);
    CSRMatrix<int> m3 = multiply_matrixCSR<int>(m1, m2);
    vector<vector<int>> multiplyResultExpected = {{0, 2, 3}, {42, 31, 12}, {63, 32, 0}};
    // check the multiply
    CHECKCSR(m3, multiplyResultExpected);
}

TEST_CASE("CSR multiply zero matrix")
{
    vector<vector<int>> array = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    vector<vector<int>> array2 = {{0, 1}, {0, -5}, {7, 8}, {56, 76}};

    CSRMatrix<int> m1 = from_vector_CSR<int>(array);
    CSRMatrix<int> m2 = from_vector_CSR<int>(array2);
    CSRMatrix<int> m3 = multiply_matrixCSR<int>(m1, m2);
    vector<vector<int>> multiplyResultExpected = {{0, 0}, {0, 0}, {0, 0}};
    // check the multiply
    CHECKCSR(m3, multiplyResultExpected);
}

TEST_CASE("CSR multiply zero row first matrix")
{
    vector<vector<int>> array = {{0, 5, 0, 3}, {0, 0, 0, 0}, {0, 7, 90, 0}};
    vector<vector<int>> array2 = {{0, 1}, {0, -5}, {7, 8}, {56, 76}};

    CSRMatrix<int> m1 = from_vector_CSR<int>(array);
    CSRMatrix<int> m2 = from_vector_CSR<int>(array2);
    CSRMatrix<int> m3 = multiply_matrixCSR<int>(m1, m2);
    vector<vector<int>> multiplyResultExpected = {{168, 203}, {0, 0}, {630, 685}};
    // check the multiply
    CHECKCSR(m3, multiplyResultExpected);
}

TEST_CASE("CSR multiply zero row second matrix")
{
    vector<vector<int>> array = {{0, 5, 0, 3}, {0, 6, 7, 0}, {0, 7, 90, 0}};
    vector<vector<int>> array2 = {{0, 0}, {0, -5}, {0, 0}, {56, 76}};

    CSRMatrix<int> m1 = from_vector_CSR<int>(array);
    CSRMatrix<int> m2 = from_vector_CSR<int>(array2);
    CSRMatrix<int> m3 = multiply_matrixCSR<int>(m1, m2);
    vector<vector<int>> multiplyResultExpected = {{168, 203}, {0, -30}, {0, -35}};
    // check the multiply
    CHECKCSR(m3, multiplyResultExpected);
}

TEST_CASE("CSR multiply zero row both matrix")
{
    vector<vector<int>> array = {{0, 5, 0, 3}, {0, 0, 0, 0}, {0, 7, 90, 0}};
    vector<vector<int>> array2 = {{0, 0}, {0, -5}, {0, 0}, {56, 76}};

    CSRMatrix<int> m1 = from_vector_CSR<int>(array);
    CSRMatrix<int> m2 = from_vector_CSR<int>(array2);
    CSRMatrix<int> m3 = multiply_matrixCSR<int>(m1, m2);
    vector<vector<int>> multiplyResultExpected = {{168, 203}, {0, 0}, {0, -35}};
    // check the multiply
    CHECKCSR(m3, multiplyResultExpected);
}

TEST_CASE("CSR multiply Exceptions")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2}, {0, -5}};

    CSRMatrix<int> m1 = from_vector_CSR<int>(array);
    CSRMatrix<int> m2 = from_vector_CSR<int>(array2);

    CHECK_THROWS_WITH_AS(multiply_matrixCSR<int>(m1, m2), "The number of columns in the first matrix must match the number of rows in the second matrix.", std::exception);
}

TEST_CASE("testing CSR subtraction")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}};
    CSRMatrix<int> m1 = from_vector_CSR<int>(array);
    CSRMatrix<int> m2 = from_vector_CSR<int>(array2);
    CSRMatrix<int> m3 = subtract_matrixCSR<int>(m1, m2);
    vector<vector<int>> subtractResultExpected = {{1, -2, -3}, {4, 10, 6}, {-7, 0, 9}};
    // check the subtract
    CHECKCSR(m3, subtractResultExpected);
}

TEST_CASE("testing CSR subtraction exceptions")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}};
    CSRMatrix<int> m1 = from_vector_CSR<int>(array);
    CSRMatrix<int> m2 = from_vector_CSR<int>(array2);
    CHECK_THROWS_WITH_AS(subtract_matrixCSR<int>(m1, m2), "The number of rows in the first matrix must match the number of rows in the second matrix.", std::exception);
}

TEST_CASE("testing CSR scalar multiply")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    CSRMatrix<int> m1 = from_vector_CSR<int>(array);
    CSRMatrix<int> m3 = scalar_multiply_CSR<int>(m1, 2);
    vector<vector<int>> scalarMultiplyResultExpected = {{2, 0, 0}, {8, 10, 12}, {0, 16, 18}};
    // check the scalar multiply
    CHECKCSR(m3, scalarMultiplyResultExpected);
}

TEST_CASE("testing CSR find max value 1")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    CSRMatrix<int> m1 = from_vector_CSR<int>(array);
    int max = find_max_CSR<int>(m1);
    CHECK(max == 9);
}

TEST_CASE("testing CSR find max value 2")
{
    vector<vector<int>> array = {{-1, 0, 0}, {-4, -5, -6}, {0, -8, -9}};
    CSRMatrix<int> m1 = from_vector_CSR<int>(array);
    int max = find_max_CSR<int>(m1);
    CHECK(max == -1);
}

TEST_CASE("testing CSR find min value 1")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    CSRMatrix<int> m1 = from_vector_CSR<int>(array);
    int min = find_min_CSR<int>(m1);
    CHECK(min == 1);
}

TEST_CASE("testing CSR find min value 2")
{
    vector<vector<int>> array = {{-1, 0, 0}, {-4, -5, -6}, {0, -8, -9}};
    CSRMatrix<int> m1 = from_vector_CSR<int>(array);
    int min = find_min_CSR<int>(m1);
    CHECK(min == -9);
}


/// @brief Ensure compressed sparse column matrix is the same as expected dense matrix
/// @param mResult The CSC matrix to compare to
/// @param mCheck The dense matrix to compare to
void CHECKCSC(CSCMatrix<int> mResult, vector<vector<int>> mCheck)
{
    CHECK(mResult.numRows == mCheck.size());
    for (size_t j = 0; j < mResult.numColumns; j++)
    {
        CHECK(mResult.numRows == mCheck[j].size());
        for (size_t i = 0; i < mResult.numRows; i++)
        {
            CHECK_MESSAGE(get_matrixCSC(mResult, i, j) == mCheck[i][j], "i = " << i << ", j = " << j);
        }
    }
}


// write similar tests for CSC
TEST_CASE("testing CSC addition")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}};
    CSCMatrix<int> m1 = from_vector_CSC<int>(array);
    CSCMatrix<int> m2 = from_vector_CSC<int>(array2);
    CSCMatrix<int> m3 = add_matrixCSC<int>(m1, m2);
    vector<vector<int>> addResultExpected = {{1, 2, 3}, {4, 0, 6}, {7, 16, 9}};
    // check the addition
    CHECKCSC(m3, addResultExpected);
}

TEST_CASE("testing CSC addition exceptions")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}, {1, 2, 3}};
    CSCMatrix<int> m1 = from_vector_CSC<int>(array);
    CSCMatrix<int> m2 = from_vector_CSC<int>(array2);
    CHECK_THROWS_WITH_AS(add_matrixCSC<int>(m1, m2), "The number of rows in the first matrix must match the number of rows in the second matrix.", std::exception);
}

TEST_CASE("testing CSC subtraction")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}};
    CSCMatrix<int> m1 = from_vector_CSC<int>(array);
    CSCMatrix<int> m2 = from_vector_CSC<int>(array2);
    CSCMatrix<int> m3 = subtract_matrixCSC<int>(m1, m2);
    vector<vector<int>> subtractResultExpected = {{1, -2, -3}, {4, 10, 6}, {-7, 0, 9}};
    // check the subtract
    CHECKCSC(m3, subtractResultExpected);
}

TEST_CASE("testing CSC subtraction exceptions")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}, {1, 2, 3}};
    CSCMatrix<int> m1 = from_vector_CSC<int>(array);
    CSCMatrix<int> m2 = from_vector_CSC<int>(array2);
    CHECK_THROWS_WITH_AS(subtract_matrixCSC<int>(m1, m2),"The number of rows in the first matrix must match the number of rows in the second matrix.", std::exception);
}

TEST_CASE("testing CSC multiplication")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}};
    CSCMatrix<int> m1 = from_vector_CSC<int>(array);
    CSCMatrix<int> m2 = from_vector_CSC<int>(array2);
    CSCMatrix<int> m3 = multiply_matrixCSC<int>(m1, m2);
    vector<vector<int>> multiplyResultExpected = {{0, 2, 3}, {42, 31, 12}, {63, 32, 0}};
    // check the multiply
    CHECKCSC(m3, multiplyResultExpected);
}

TEST_CASE("testing CSC multiplication exceptions")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}, {1, 2, 3}};
    CSCMatrix<int> m1 = from_vector_CSC<int>(array);
    CSCMatrix<int> m2 = from_vector_CSC<int>(array2);
    CHECK_THROWS_WITH_AS(multiply_matrixCSC<int>(m1, m2), "The number of columns in the first matrix must match the number of rows in the second matrix.", std::exception);
}

TEST_CASE("testing CSC find max value 1")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    CSCMatrix<int> m1 = from_vector_CSC<int>(array);
    int max = find_max_CSC<int>(m1);
    CHECK(max == 9);
}

TEST_CASE("testing CSC find max value 2")
{
    vector<vector<int>> array = {{-1, 0, 0}, {-4, -5, -6}, {0, -8, -9}};
    CSCMatrix<int> m1 = from_vector_CSC<int>(array);
    int max = find_max_CSC<int>(m1);
    CHECK(max == -1);
}

TEST_CASE("testing CSC find min value 1")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    CSCMatrix<int> m1 = from_vector_CSC<int>(array);
    int min = find_min_CSC<int>(m1);
    CHECK(min == 1);
}

TEST_CASE("testing CSC find min value 2")
{
    vector<vector<int>> array = {{-1, 0, 0}, {-4, -5, -6}, {0, -8, -9}};
    CSCMatrix<int> m1 = from_vector_CSC<int>(array);
    int min = find_min_CSC<int>(m1);
    CHECK(min == -9);
}


TEST_CASE("testing CSC scalar multiplication")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    CSCMatrix<int> m1 = from_vector_CSC<int>(array);
    CSCMatrix<int> m2 = scalar_multiply_CSC<int>(m1, 2);
    vector<vector<int>> scalarMultiplyResultExpected = {{2, 0, 0}, {8, 10, 12}, {0, 16, 18}};
    // check the scalar multiply
    CHECKCSC(m2, scalarMultiplyResultExpected);
}

// test for Gaussian Elimination
TEST_CASE("testing Gaussian Elimination 1")
{
    std::vector<std::vector<double>> A1 = {{2, 1, -1},
                                           {-3, -1, 2},
                                           {-2, 1, 2}};
    std::vector<double> b1 = {8, -11, -3};
    std::vector<double> x1 = {2, 3, -1};

    bool result1 = gaussian_elimination(A1, b1);
    CHECK(result1);
    CHECK_VECTOR_EQ(b1, x1, 1e-6);
}

TEST_CASE("testing Gaussian Elimination 2")
{
    std::vector<std::vector<double>> A2 = {{1, 2, -1},
                                           {2, 1, -2},
                                           {-3, 1, 1}};
    std::vector<double> b2 = {3, 3, -6};
    std::vector<double> x2 = {3, 1, 2};

    bool result2 = gaussian_elimination(A2, b2);
    CHECK(result2);
    CHECK_VECTOR_EQ(b2, x2, 1e-6);
}

TEST_CASE("testing Gaussian Elimination 3")
{
    std::vector<std::vector<double>> A3 = {{1, 2},
                                           {2, 4}};
    std::vector<double> b3 = {3, 6};
    bool result3 = gaussian_elimination(A3, b3);
    CHECK(!result3);
}

// test for QR decomposition
TEST_CASE("QR Factorization Test 1")
{
    // Create a matrix
    std::vector<std::vector<double>> A = {{1.0, 2.0, 3.0},
                                          {4.0, 5.0, 6.0},
                                          {7.0, 8.0, 7.0},
                                          {4.0, 2.0, 1.0}};

    // Call the QR factorization function
    auto qr = qr_factorization(A);
    auto Q = qr.first;
    auto R = qr.second;

    // Check that the dimensions of the Q and R matrices are correct
    CHECK(Q.size() == A.size());
    CHECK(Q[0].size() == A[0].size());
    CHECK(R.size() == A[0].size());
    CHECK(R[0].size() == A[0].size());

    // Check that Q is orthogonal
    std::vector<std::vector<double>> QT = transpose(Q);
    std::vector<std::vector<double>> Q_QT = mult_matrix(QT, Q);
    std::vector<std::vector<double>> I = identity_matrix(A[0].size());
    
    CHECK_MATRIX_EQ(Q_QT, I, 1e-6);

    // Check that Q*R = A
    std::vector<std::vector<double>> Q_times_R = mult_matrix(Q, R);
    for (size_t i = 0; i < A.size(); i++)
    {
        for (size_t j = 0; j < A[0].size(); j++)
        {
            CHECK(doctest::Approx(Q_times_R[i][j]).epsilon(1e-6) == A[i][j]);
        }
    }
}

// test for LU factorization
TEST_CASE("LU Factorization test 1")
{
    // Test a 2x2 matrix
    std::vector<std::vector<double>> A = {{4, 3}, {6, 3}};

    auto result = lu_factorization(A);
    std::vector<std::vector<double>> P = std::get<0>(result);
    std::vector<std::vector<double>> L = std::get<1>(result);
    std::vector<std::vector<double>> U = std::get<2>(result);

    // Check that L*U = P*A
    std::vector<std::vector<double>> L_U = mult_matrix(L, U);
    std::vector<std::vector<double>> P_A = mult_matrix(P, A);
    CHECK_MATRIX_EQ(L_U, P_A, 1e-6);
}

TEST_CASE("LU Factorization test 2")
{
    std::vector<std::vector<double>> A = {{1, 2, -1, 4}, {-2, -3, 4, 5}, 
                                          {3, 6, -2, 7}, {1, 3, 1, 9}};
    
    auto result = lu_factorization(A);
    std::vector<std::vector<double>> P = std::get<0>(result);
    std::vector<std::vector<double>> L = std::get<1>(result);
    std::vector<std::vector<double>> U = std::get<2>(result);

    // Check that L*U = P*A
    std::vector<std::vector<double>> L_U = mult_matrix(L, U);
    std::vector<std::vector<double>> P_A = mult_matrix(P, A);
    CHECK_MATRIX_EQ(L_U, P_A, 1e-6);
}

// test for LU factorization
TEST_CASE("LU Factorization in Place")
{
    // Test a 2x2 matrix
    std::vector<std::vector<double>> A = {{1, 2}, {3, 4}};

    auto result = lu_factorization_inplace(A);

    std::vector<std::vector<double>> expected = {{3, 4}, {0.3333333, 0.66666667}};
        CHECK_MATRIX_EQ(A, expected, 1e-6);
}


// test for LDL factorization
TEST_CASE("LDL^T Factorization")
{
    std::vector<std::vector<double>> A = {{4, 12, -16}, {12, 37, -43}, {-16, -43, 98}};
    auto result = ldlt_factorization(A);
    std::vector<std::vector<double>> L = result.first;
    std::vector<double> d = result.second;
    std::vector<std::vector<double>> D(A.size(), std::vector<double>(A.size(), 0));
    for (size_t i = 0; i < A.size(); i++) {D[i][i] = d[i];}

    // Check L is correct
    std::vector<std::vector<double>> L_correct = {{1, 0, 0}, {3, 1, 0}, {-4, 5, 1}};
    CHECK_MATRIX_EQ(L, L_correct, 1e-6);

    // Check D is correct
    std::vector<double> d_correct = {{4, 1, 9}};
    CHECK_VECTOR_EQ(d, d_correct, 1e-6);

    // Check that LDL^T = A
    std::vector<std::vector<double>> LD = mult_matrix(L, D);
    std::vector<std::vector<double>> LT = transpose(L);
    std::vector<std::vector<double>> LDLT = mult_matrix(LD, LT);
    CHECK_MATRIX_EQ(LDLT, A, 1e-6);
}

// test for Cholesky factorization
TEST_CASE("Cholesky Factorization 0")
{
    std::vector<std::vector<double>> A = {{2, -1, 0},
                                          {-1, 2, -1},
                                          {0, -1, 2}};
    std::vector<std::vector<double>> L = cholesky_factorization(A);
    //PRINT_MATRIX(L);

    // Check that L*L^T = A
    std::vector<std::vector<double>> LT = transpose(L);
    std::vector<std::vector<double>> LLT = mult_matrix(L, LT);
    CHECK_MATRIX_EQ(LLT, A, 1e-6);
}


TEST_CASE("Cholesky Factorization 1")
{
    std::vector<std::vector<double>> A = {{4, 12, -16}, {12, 37, -43}, {-16, -43, 98}};
    std::vector<std::vector<double>> L = cholesky_factorization(A);

    // Check that L*L^T = A
    std::vector<std::vector<double>> LT = transpose(L);
    std::vector<std::vector<double>> LLT = mult_matrix(L, LT);
    CHECK_MATRIX_EQ(LLT, A, 1e-6);
}

TEST_CASE("Gauss-Seidel")
{
    std::vector<std::vector<double>> A = {{16, 3},
                                          {7, -11}};
    std::vector<double> b = {19, -4};
    std::vector<double> x = {0, 0};

    const double tol = 1e-6;
    int max_iter = 100;

    std::vector<double> x_expected = {1.0, 1.0};
    bool success = gauss_seidel(A, b, x, max_iter);

    CHECK(success);
    CHECK_VECTOR_EQ(x, x_expected, tol);
}

TEST_CASE("Jacobi Iteration 0")
{
    const std::vector<std::vector<double>> A = {{4.0, 1.0, 1.0}, {1.0, 4.0, 1.0}, {1.0, 1.0, 4.0}};
    const std::vector<double> b = {6.0, 6.0, 6.0};
    const double tol = 1e-6;
    const int max_iter = 100;

    std::vector<double> x_expected = {1.0, 1.0, 1.0};
    std::vector<double> x = jacobi_iteration(A, b, tol, max_iter);

    CHECK_VECTOR_EQ(x, x_expected, tol);
}

TEST_CASE("Jacobi Iteration 1")
{
    const std::vector<std::vector<double>> A = {{3.0, 1.0, 1.0}, {1.0, 5.0, 2.0}, {2.0, 3.0, 6.0}};
    const std::vector<double> b = {5.0, 10.0, 15.0};
    const double tol = 1e-6;
    const int max_iter = 100;

    std::vector<double> x_expected = {0.714286, 1.190476, 1.666667};
    std::vector<double> x = jacobi_iteration(A, b, tol, max_iter);

    CHECK_VECTOR_EQ(x, x_expected, tol);
}


TEST_CASE("Jacobi Iteration CSR 0")
{
    std::vector<std::vector<double>> A = {{4.0, 1.0, 1.0}, {1.0, 4.0, 1.0}, {1.0, 1.0, 4.0}};
    auto CSR_A = from_vector_CSR<double>(A);
    std::vector<double> b = {6.0, 6.0, 6.0};
    const double tol = 1e-6;
    const int max_iter = 100;

    std::vector<double> x_expected = {1.0, 1.0, 1.0};
    std::vector<double> x = jacobi_method_CSR<double>(CSR_A, b,tol, max_iter);

    CHECK_VECTOR_EQ(x, x_expected, tol);
}

TEST_CASE("Jacobi Iteration CSR 1")
{
    std::vector<std::vector<double>> A = {{3.0, 1.0, 1.0}, {1.0, 5.0, 2.0}, {2.0, 3.0, 6.0}};
    auto CSR_A = from_vector_CSR<double>(A);
    const std::vector<double> b = {5.0, 10.0, 15.0};
    const double tol = 1e-6;
    const int max_iter = 100;

    std::vector<double> x_expected = {0.714286, 1.190476, 1.666667};
    std::vector<double> x = jacobi_method_CSR(CSR_A, b, tol, max_iter);

    CHECK_VECTOR_EQ(x, x_expected, tol);
}

TEST_CASE("SSOR Iteration (w = 1.0)")
{
    const std::vector<std::vector<double>> A = {{4.0, 1.0, 1.0}, {1.0, 4.0, 1.0}, {1.0, 1.0, 4.0}};
    const std::vector<double> b = {6.0, 6.0, 6.0};
    const double tol = 1e-6;
    const int max_iter = 100;
    const double omega = 1.0;  // Set the relaxation parameter (w) to 1.0

    std::vector<double> x_expected = {1.0, 1.0, 1.0};
    std::vector<double> x = ssor_iteration(A, b, tol, max_iter, omega);

    CHECK_VECTOR_EQ(x, x_expected, tol);
}

TEST_CASE("SSOR Iteration (w = 1.5)")
{
    const std::vector<std::vector<double>> A = {{4.0, 1.0, 1.0}, {1.0, 4.0, 1.0}, {1.0, 1.0, 4.0}};
    const std::vector<double> b = {6.0, 6.0, 6.0};
    const double tol = 1e-6;
    const int max_iter = 100;
    const double omega = 1.5;  // Set the relaxation parameter (w) to 1.5

    std::vector<double> x_expected = {1.0, 1.0, 1.0};
    std::vector<double> x = ssor_iteration(A, b, tol, max_iter, omega);

    CHECK_VECTOR_EQ(x, x_expected, tol);
}

TEST_CASE("SSOR Iteration (w = 0.5)")
{
    const std::vector<std::vector<double>> A = {{4.0, 1.0, 1.0}, {1.0, 4.0, 1.0}, {1.0, 1.0, 4.0}};
    const std::vector<double> b = {6.0, 6.0, 6.0};
    const double tol = 1e-2;
    const int max_iter = 100;
    const double omega = 0.5;  // Set the relaxation parameter (w) to 0.5

    std::vector<double> x_expected = {1.0, 1.0, 1.0};
    std::vector<double> x = ssor_iteration(A, b, tol, max_iter, omega);

    CHECK_VECTOR_EQ(x, x_expected, tol);
}

// TEST_CASE("Incomplete Cholesky Factorization 1") {
//     // Create a matrix A
//     vector<vector<double>> A = {{2, -1, 0},
//                                 {-1, 2, -1},
//                                 {0, -1, 2}};

//     // Compute the incomplete Cholesky factorization of A
//     vector<vector<double>> K = incompleteCholesky(A, 1e-1);
//     // PRINT_MATRIX(K);
//     auto KKt = mult_matrix(K, transpose(K));
//     CHECK_MATRIX_EQ(KKt, A, 1e-12);
// }

// TEST_CASE("Incomplete Cholesky Factorization 2") {
//     // Create a matrix A
//     vector<vector<double>> A = {{4, 0, 0, 0}, 
//                                 {0, 6, 0, 2}, 
//                                 {0, 0, 8, 0}, 
//                                 {0, 2, 0, 10}};

//     // Compute the incomplete Cholesky factorization of A
//     vector<vector<double>> K = incompleteCholesky(A, 1e-12);
//     auto KKt = mult_matrix(K, transpose(K));
//     vector<vector<double>> KKt_expected = {{4, 0, 0, 0}, 
//                                            {0, 6, 0, 2}, 
//                                            {0, 0, 8, 0}, 
//                                            {0, 2, 0, 10}};
//     CHECK_MATRIX_EQ(KKt, A, 1e-12);
// }

// TEST_CASE("Incomplete Cholesky Factorization 3") {
//     std::vector<std::vector<double>> A = {{4, 12, -16}, {12, 37, -43}, {-16, -43, 98}};
//     std::vector<std::vector<double>> K = incompleteCholesky(A, 1e-12);
//     auto KKt = mult_matrix(K, transpose(K));
//     CHECK_MATRIX_EQ(KKt, A, 1e-12);
// }

// TEST_CASE("Incomplete Cholesky Factorization 4") {
//     // Create a matrix A
//     vector<vector<double>> A = {{3.0, -1.0, 0.0},
//                                 {-1.0, 3.0, -1.0},
//                                 {0.0, -1.0, 3.0}};

//     // Compute the incomplete Cholesky factorization of A
//     vector<vector<double>> K = incompleteCholesky(A, 1e-12);
//     // PRINT_MATRIX(K);
//     auto KKt = mult_matrix(K, transpose(K));
//     CHECK_MATRIX_EQ(KKt, A, 1e-12);
// }

// TEST_CASE("Matrix Inverse 1") {
//     // Create a matrix A
//     vector<vector<double>> A = {{6.0,2.0,3.0}, 
//                                 {1.0,1.0,1.0}, 
//                                 {0.0,4.0,9.0}};

//     // Compute the incomplete Cholesky factorization of A
//     vector<vector<double>> A_inverse = {{0.208333,-0.25,-0.0416667}, 
//                                         {-0.375,  2.25, -0.125}, 
//                                         {0.166667, -1.0, 0.166667}};
//     // PRINT_MATRIX(K);
//     vector<vector<double>> test_inverse = matrix_inverse(A);
    
//     CHECK_MATRIX_EQ(A_inverse, test_inverse, 1e-6);
// }

TEST_CASE("dense ADD load file corectness") {
    std::vector<std::vector<double>> m1 = load_fileMatrix<double>("../../../data/matrices/small_test_matrix.mtx");
    std::vector<std::vector<double>> m2 = load_fileMatrix<double>("../../../data/matrices/small_test_matrix.mtx");
    std::vector<std::vector<double>> m3 = sum_matrix(m1, m2);
    vector<vector<double>> expected = {{1.0, 0.0, 4.0},
                                {2.0, 0.0, -8.0},
                                {6.0, 0.0, 0.0}};

    CHECKMATRIX(m3,expected);
}

TEST_CASE("CSR ADD load file corectness") {
    CSRMatrix<double> m1 = load_fileCSR<double>("../../../data/matrices/small_test_matrix.mtx");
    CSRMatrix<double> m2 = load_fileCSR<double>("../../../data/matrices/small_test_matrix.mtx");

    CSRMatrix<double> m3 = add_matrixCSR<double>(m1, m2);
    vector<vector<double>> expected = {{1.0, 0.0, 4.0},
                                {2.0, 0.0, -8.0},
                                {6.0, 0.0, 0.0}};
    CHECKCSR(m3,expected);
}

TEST_CASE("dense ADD load file corectness two") {
    std::vector<std::vector<double>> m1 = load_fileMatrix<double>("../../../data/matrices/small_test_matrix.mtx");
    std::vector<std::vector<double>> m2 = load_fileMatrix<double>("../../../data/matrices/small_test_matrix_two.mtx");
    vector<vector<double>> expectedOne = {{0.5, 0.0, 2.0},
                                {1.0, 0.0, -4.0},
                                {3.0, 0.0, 0.0}};
    CHECKMATRIX(m1,expectedOne);
    vector<vector<double>> expectedTwo = {{1.0, 1.0, 1.0},
                                {1.0, 1.0, 1.0},
                                {1.0, 1.0, 1.0}};
    CHECKMATRIX(m2,expectedTwo);
    std::vector<std::vector<double>> m3 = sum_matrix(m1, m2);
    vector<vector<double>> expectedThree = {{1.5, 1.0, 3.0},
                                {2.0, 1.0, -3.0},
                                {4.0, 1.0, 1.0}};
    CHECKMATRIX(m3,expectedThree);
    std::vector<std::vector<double>> m4 = sub_matrix(m1, m2);
    vector<vector<double>> expectedFour = {{-0.5, -1.0, 1.0},
                                {0.0, -1.0, -5.0},
                                {2.0, -1.0, -1.0}};
    CHECKMATRIX(m4,expectedFour);
}

TEST_CASE("CSR ADD load file corectness two") {
    CSRMatrix<double> m1 = load_fileCSR<double>("../../../data/matrices/small_test_matrix.mtx");
    CSRMatrix<double> m2 = load_fileCSR<double>("../../../data/matrices/small_test_matrix_two.mtx");
    vector<vector<double>> expectedOne = {{0.5, 0.0, 2.0},
                                {1.0, 0.0, -4.0},
                                {3.0, 0.0, 0.0}};
    CHECKCSR(m1,expectedOne);
    vector<vector<double>> expectedTwo = {{1.0, 1.0, 1.0},
                                {1.0, 1.0, 1.0},
                                {1.0, 1.0, 1.0}};
    CHECKCSR(m2,expectedTwo);
    CSRMatrix<double> m3 = add_matrixCSR(m1, m2);
    vector<vector<double>> expectedThree = {{1.5, 1.0, 3.0},
                                {2.0, 1.0, -3.0},
                                {4.0, 1.0, 1.0}};
    CHECKCSR(m3,expectedThree);
    CSRMatrix<double> m4 = subtract_matrixCSR(m1, m2);
    vector<vector<double>> expectedFour = {{-0.5, -1.0, 1.0},
                                {0.0, -1.0, -5.0},
                                {2.0, -1.0, -1.0}};
    CHECKCSR(m4,expectedFour);
}

TEST_CASE("Gauss Seidel Iteration CSR 0")
{
    std::vector<std::vector<double>> A = {{4.0, 1.0, 1.0}, {1.0, 4.0, 1.0}, {1.0, 1.0, 4.0}};
    auto CSR_A = from_vector_CSR<double>(A);
    std::vector<double> b = {6.0, 6.0, 6.0};
    const double tol = 1e-6;
    const int max_iter = 100;

    std::vector<double> x_expected = {1.0, 1.0, 1.0};
    std::vector<double> x = gauss_sidel_CSR<double>(CSR_A, b,tol, max_iter);

    CHECK_VECTOR_EQ(x, x_expected, tol);
}

TEST_CASE("Gauss Seidel Iteration CSR 1")
{
    std::vector<std::vector<double>> A = {{3.0, 1.0, 1.0}, {1.0, 5.0, 2.0}, {2.0, 3.0, 6.0}};
    auto CSR_A = from_vector_CSR<double>(A);
    const std::vector<double> b = {5.0, 10.0, 15.0};
    const double tol = 1e-6;
    const int max_iter = 100;

    std::vector<double> x_expected = {0.714286, 1.190476, 1.666667};
    std::vector<double> x = gauss_sidel_CSR(CSR_A, b, tol, max_iter);

    CHECK_VECTOR_EQ(x, x_expected, tol);
}

TEST_CASE("SSOR Iteration 1")
{

    const std::vector<std::vector<double>> A = {{4.0, 1.0, 1.0}, {1.0, 4.0, 1.0}, {1.0, 1.0, 4.0}};
    const std::vector<double> b = {6.0, 6.0, 6.0};
    const double tol = 1e-4;
    const int max_iter = 100;
    std::vector<double> x_expected = {1.0, 1.0, 1.0};

    double omega = 0.5;  // Set the relaxation parameter (w) to 0.5
    std::vector<double> x = ssor_iteration(A, b, tol, max_iter, omega);
    CHECK_VECTOR_EQ(x, x_expected, tol);

    omega = 1.0;  // Set the relaxation parameter (w) to 1.0
    x = ssor_iteration(A, b, tol, max_iter, omega);
    CHECK_VECTOR_EQ(x, x_expected, tol);

    omega = 1.5;  // Set the relaxation parameter (w) to 1.5
    x = ssor_iteration(A, b, tol, max_iter, omega);
    CHECK_VECTOR_EQ(x, x_expected, tol);
}

TEST_CASE("SSOR Iteration 2")
{
    std::vector<std::vector<double>> A = {{16, 3},
                                          {7, 11}};
    std::vector<double> b = {19, 18};
    const double tol = 1e-2;
    const int max_iter = 100;
    std::vector<double> x_expected = {1.0, 1.0};

    double omega = 0.5;  // Set the relaxation parameter (w) to 0.5
    std::vector<double> x = ssor_iteration(A, b, tol, max_iter, omega);
    CHECK_VECTOR_EQ(x, x_expected, tol);

    omega = 1.0;  // Set the relaxation parameter (w) to 1.0
    x = ssor_iteration(A, b, tol, max_iter, omega);
    CHECK_VECTOR_EQ(x, x_expected, tol);

    omega = 1.5;  // Set the relaxation parameter (w) to 1.5
    x = ssor_iteration(A, b, tol, max_iter, omega);
    CHECK_VECTOR_EQ(x, x_expected, tol);
}

TEST_CASE("CSR SSOR Iteration 1")
{

    std::vector<std::vector<double>> arr = {{4.0, 1.0, 1.0}, {1.0, 4.0, 1.0}, {1.0, 1.0, 4.0}};
    CSRMatrix<double> A = from_vector_CSR<double>(arr);
    std::vector<double> b = {6.0, 6.0, 6.0};
    const double tol = 1e-4;
    const int max_iter = 100;
    std::vector<double> x_expected = {1.0, 1.0, 1.0};

    double omega = 0.5;  // Set the relaxation parameter (w) to 0.5
    std::vector<double> x = ssor_iteration_CSR(A, b, tol, max_iter, omega);
    CHECK_VECTOR_EQ(x, x_expected, tol);

    omega = 1.0;  // Set the relaxation parameter (w) to 1.0
    x = ssor_iteration_CSR(A, b, tol, max_iter, omega);
    CHECK_VECTOR_EQ(x, x_expected, tol);

    omega = 1.5;  // Set the relaxation parameter (w) to 1.5
    x = ssor_iteration_CSR(A, b, tol, max_iter, omega);
    CHECK_VECTOR_EQ(x, x_expected, tol);
}

TEST_CASE("CSR SSOR Iteration 2")
{

    std::vector<std::vector<double>> arr = {{16, 3},
                                          {7, 11}};
    CSRMatrix<double> A = from_vector_CSR<double>(arr);
    std::vector<double> b = {19, 18};
    const double tol = 1e-2;
    const int max_iter = 100;
    std::vector<double> x_expected = {1.0, 1.0};

    double omega = 0.5;  // Set the relaxation parameter (w) to 0.5
    std::vector<double> x = ssor_iteration_CSR(A, b, tol, max_iter, omega);
    CHECK_VECTOR_EQ(x, x_expected, tol);

    omega = 1.0;  // Set the relaxation parameter (w) to 1.0
    x = ssor_iteration_CSR(A, b, tol, max_iter, omega);
    CHECK_VECTOR_EQ(x, x_expected, tol);

    omega = 1.5;  // Set the relaxation parameter (w) to 1.5
    x = ssor_iteration_CSR(A, b, tol, max_iter, omega);
    CHECK_VECTOR_EQ(x, x_expected, tol);
}

TEST_CASE("ILU Factorization 1")
{
    std::vector<std::vector<double>> A = {{5, -2, -2, 0},
                                          {-2, 5, 0, -2},
                                          {-2, 0, 5, -2},
                                          {0, -2, -2, 5}};
    ilu(A, 0);
    // Taken from https://www.mathworks.com/help/matlab/ref/ilu.html
    std::vector<std::vector<double>> LU_expected = {{5.0000, -2.0000, -2.0000, 0.0000},
                                                    {-0.4000, 4.2000, 0.0000, -2.000},
                                                    {-0.4000, -0.0000, 4.2000, -2.000},
                                                    {0.0000, -0.4762, -0.4762, 3.0952}};
    CHECK_MATRIX_EQ(A, LU_expected, 1e-4);

}

TEST_CASE("ILUT Factorization 1")
{
    std::vector<std::vector<double>> A = {{5, -2, -2, 0},
                                          {-2, 5, 0, -2},
                                          {-2, 0, 5, -2},
                                          {0, -2, -2, 5}};
    ilut(A, 0.01);
    // Taken from https://www.mathworks.com/help/matlab/ref/ilu.html
    std::vector<std::vector<double>> LU_expected = {{5.0000, -2.0000, -2.0000, 0.0000},
                                                    {-0.4000, 4.2000, -0.8000, -2.0000},
                                                    {-0.4000, -0.1905, 4.0476, -2.3810},
                                                    {0.0000, -0.4762, -0.5882, 2.6471}};
    CHECK_MATRIX_EQ(A, LU_expected, 1e-4);

}

TEST_CASE("ILUT Factorization 2")
{
    std::vector<std::vector<double>> A = {{5, -2, -2, 0},
                                          {-2, 5, 0, -2},
                                          {-2, 0, 5, -2},
                                          {0, -2, -2, 5}};
    ilut(A, 0.035);
    std::vector<std::vector<double>> LU_expected = {{5.0000, -2.0000, -2.0000, 0.0000},
                                                    {-0.4000, 4.2000, -0.8000, -2.000},
                                                    {-0.4000, 0.0000, 4.2000, -2.0000},
                                                    {0.0000, -0.4762, -0.5669, 2.9138}};
    CHECK_MATRIX_EQ(A, LU_expected, 1e-4);
}

TEST_CASE("GCR Iteration 1")
{
    const std::vector<std::vector<double>> A = {{4.0, 1.0, 1.0}, {1.0, 4.0, 1.0}, {1.0, 1.0, 4.0}};
    const std::vector<double> b = {6.0, 6.0, 6.0};
    const double tol = 1e-4;
    const int max_iter = 100;
    std::vector<double> x_initial = {0.0, 0.0, 0.0};
    std::vector<double> x_expected = {1.0, 1.0, 1.0};

    std::vector<double> x = gcr(A, b, x_initial, tol, max_iter);
    CHECK_VECTOR_EQ(x, x_expected, tol);
}