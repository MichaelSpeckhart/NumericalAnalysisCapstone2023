#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../functionsCSRParallel.cc"
#include "../functions.cc"
#include "../functionsCSR.cc"
#include "../functionsCSC.cc"
#include "fstream"
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

TEST_CASE("testing CSR Add")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};

    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    // check m1
    CHECKCSR(m1, array);

    CSRMatrix<int> m2 = from_vector<int>(array2);
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

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);

    CHECK_THROWS_WITH_AS(add_matrixCSR<int>(m1, m2), "The number of rows in the first matrix must match the number of rows in the second matrix.", std::exception);

    vector<vector<int>> array3 = {{1, 0}, {4, 5}, {0, 8}};
    CSRMatrix<int> m3 = from_vector<int>(array3);

    CHECK_THROWS_WITH_AS(add_matrixCSR<int>(m1, m3), "The number of columns in the first matrix must match the number of columns in the second matrix.", std::exception);
}

TEST_CASE("testing CSR transpose")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};

    CSRMatrix<int> m1 = from_vector<int>(array);

    CSRMatrix<int> m2 = transpose_matrixCSR<int>(m1);
    vector<vector<int>> transposeResultExpected = {{1, 4, 0}, {0, 5, 8}, {0, 6, 9}};
    // check the transpose
    CHECKCSR(m2, transposeResultExpected);
}

TEST_CASE("testing CSR multiply")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = multiply_matrixCSR<int>(m1, m2);
    vector<vector<int>> multiplyResultExpected = {{0, 2, 3}, {42, 31, 12}, {63, 32, 0}};
    // check the multiply
    CHECKCSR(m3, multiplyResultExpected);
}

TEST_CASE("CSR multiply zero matrix")
{
    vector<vector<int>> array = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    vector<vector<int>> array2 = {{0, 1}, {0, -5}, {7, 8}, {56, 76}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = multiply_matrixCSR<int>(m1, m2);
    vector<vector<int>> multiplyResultExpected = {{0, 0}, {0, 0}, {0, 0}};
    // check the multiply
    CHECKCSR(m3, multiplyResultExpected);
}

TEST_CASE("CSR multiply zero row first matrix")
{
    vector<vector<int>> array = {{0, 5, 0, 3}, {0, 0, 0, 0}, {0, 7, 90, 0}};
    vector<vector<int>> array2 = {{0, 1}, {0, -5}, {7, 8}, {56, 76}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = multiply_matrixCSR<int>(m1, m2);
    vector<vector<int>> multiplyResultExpected = {{168, 203}, {0, 0}, {630, 685}};
    // check the multiply
    CHECKCSR(m3, multiplyResultExpected);
}

TEST_CASE("CSR multiply zero row second matrix")
{
    vector<vector<int>> array = {{0, 5, 0, 3}, {0, 6, 7, 0}, {0, 7, 90, 0}};
    vector<vector<int>> array2 = {{0, 0}, {0, -5}, {0, 0}, {56, 76}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = multiply_matrixCSR<int>(m1, m2);
    vector<vector<int>> multiplyResultExpected = {{168, 203}, {0, -30}, {0, -35}};
    // check the multiply
    CHECKCSR(m3, multiplyResultExpected);
}

TEST_CASE("CSR multiply zero row both matrix")
{
    vector<vector<int>> array = {{0, 5, 0, 3}, {0, 0, 0, 0}, {0, 7, 90, 0}};
    vector<vector<int>> array2 = {{0, 0}, {0, -5}, {0, 0}, {56, 76}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = multiply_matrixCSR<int>(m1, m2);
    vector<vector<int>> multiplyResultExpected = {{168, 203}, {0, 0}, {0, -35}};
    // check the multiply
    CHECKCSR(m3, multiplyResultExpected);
}

TEST_CASE("CSR multiply Exceptions")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2}, {0, -5}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);

    CHECK_THROWS_WITH_AS(multiply_matrixCSR<int>(m1, m2), "The number of columns in the first matrix must match the number of rows in the second matrix.", std::exception);
}

TEST_CASE("testing CSR subtraction")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = subtract_matrixCSR<int>(m1, m2);
    vector<vector<int>> subtractResultExpected = {{1, -2, -3}, {4, 10, 6}, {-7, 0, 9}};
    // check the subtract
    CHECKCSR(m3, subtractResultExpected);
}

TEST_CASE("testing CSR subtraction exceptions")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CHECK_THROWS_WITH_AS(subtract_matrixCSR<int>(m1, m2), "The number of rows and columns in the first matrix must match the number of rows and columns in the second matrix.", std::exception);
}

TEST_CASE("testing CSR scalar multiply")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m3 = scalar_multiply_matrixCSR<int>(m1, 2);
    vector<vector<int>> scalarMultiplyResultExpected = {{2, 0, 0}, {8, 10, 12}, {0, 16, 18}};
    // check the scalar multiply
    CHECKCSR(m3, scalarMultiplyResultExpected);
}

TEST_CASE("testing CSR find max value 1")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    int max = find_max_valueCSR<int>(m1);
    CHECK(max == 9);
}

TEST_CASE("testing CSR find max value 2")
{
    vector<vector<int>> array = {{-1, 0, 0}, {-4, -5, -6}, {0, -8, -9}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    int max = find_max_valueCSR<int>(m1);
    CHECK(max == -1);
}

TEST_CASE("testing CSR find min value 1")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    int min = find_min_valueCSR<int>(m1);
    CHECK(min == 1);
}

TEST_CASE("testing CSR find min value 2")
{
    vector<vector<int>> array = {{-1, 0, 0}, {-4, -5, -6}, {0, -8, -9}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    int min = find_min_valueCSR<int>(m1);
    CHECK(min == -9);
}

TEST_CASE("testing CSR find max value exception")
{
    vector<vector<int>> array = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    CHECK_THROWS_WITH_AS(find_max_valueCSR<int>(m1), "The matrix is empty.", std::exception);
}

TEST_CASE("testing CSR find min value exception")
{
    vector<vector<int>> array = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    CHECK_THROWS_WITH_AS(find_min_valueCSR<int>(m1), "The matrix is empty.", std::exception);
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
    CSCMatrix<int> m1 = from_vector<int>(array);
    CSCMatrix<int> m2 = from_vector<int>(array2);
    CSCMatrix<int> m3 = add_matrixCSC<int>(m1, m2);
    vector<vector<int>> addResultExpected = {{1, 2, 3}, {4, 0, 6}, {7, 16, 9}};
    // check the addition
    CHECKCSC(m3, addResultExpected);
}

TEST_CASE("testing CSC addition exceptions")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}, {1, 2, 3}};
    CSCMatrix<int> m1 = from_vector<int>(array);
    CSCMatrix<int> m2 = from_vector<int>(array2);
    CHECK_THROWS_WITH_AS(add_matrixCSC<int>(m1, m2), "The number of rows and columns in the first matrix must match the number of rows and columns in the second matrix.", std::exception);
}

TEST_CASE("testing CSC subtraction")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}};
    CSCMatrix<int> m1 = from_vector<int>(array);
    CSCMatrix<int> m2 = from_vector<int>(array2);
    CSCMatrix<int> m3 = subtract_matrixCSC<int>(m1, m2);
    vector<vector<int>> subtractResultExpected = {{1, -2, -3}, {4, 10, 6}, {-7, 0, 9}};
    // check the subtract
    CHECKCSC(m3, subtractResultExpected);
}

TEST_CASE("testing CSC subtraction exceptions")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}, {1, 2, 3}};
    CSCMatrix<int> m1 = from_vector<int>(array);
    CSCMatrix<int> m2 = from_vector<int>(array2);
    CHECK_THROWS_WITH_AS(subtract_matrixCSC<int>(m1, m2), "The number of rows and columns in the first matrix must match the number of rows and columns in the second matrix.", std::exception);
}

TEST_CASE("testing CSC multiplication")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}};
    CSCMatrix<int> m1 = from_vector<int>(array);
    CSCMatrix<int> m2 = from_vector<int>(array2);
    CSCMatrix<int> m3 = multiply_matrixCSC<int>(m1, m2);
    vector<vector<int>> multiplyResultExpected = {{0, 2, 3}, {0, -10, 0}, {0, 16, 0}};
    // check the multiply
    CHECKCSC(m3, multiplyResultExpected);
}

TEST_CASE("testing CSC multiplication exceptions")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    vector<vector<int>> array2 = {{0, 2, 3}, {0, -5, 0}, {7, 8, 0}, {1, 2, 3}};
    CSCMatrix<int> m1 = from_vector<int>(array);
    CSCMatrix<int> m2 = from_vector<int>(array2);
    CHECK_THROWS_WITH_AS(multiply_matrixCSC<int>(m1, m2), "The number of columns in the first matrix must match the number of rows in the second matrix.", std::exception);
}

TEST_CASE("testing CSC find max value 1")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    CSCMatrix<int> m1 = from_vector<int>(array);
    int max = find_max_valueCSC<int>(m1);
    CHECK(max == 9);
}

TEST_CASE("testing CSC find max value 2")
{
    vector<vector<int>> array = {{-1, 0, 0}, {-4, -5, -6}, {0, -8, -9}};
    CSCMatrix<int> m1 = from_vector<int>(array);
    int max = find_max_valueCSC<int>(m1);
    CHECK(max == -1);
}

TEST_CASE("testing CSC find min value 1")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    CSCMatrix<int> m1 = from_vector<int>(array);
    int min = find_min_valueCSC<int>(m1);
    CHECK(min == 1);
}

TEST_CASE("testing CSC find min value 2")
{
    vector<vector<int>> array = {{-1, 0, 0}, {-4, -5, -6}, {0, -8, -9}};
    CSCMatrix<int> m1 = from_vector<int>(array);
    int min = find_min_valueCSC<int>(m1);
    CHECK(min == -9);
}

TEST_CASE("testing CSC find max value exceptions")
{
    vector<vector<int>> array = {{}};
    CSCMatrix<int> m1 = from_vector<int>(array);
    CHECK_THROWS_WITH_AS(find_max_valueCSC<int>(m1), "The matrix is empty.", std::exception);
}

TEST_CASE("testing CSC find min value exceptions")
{
    vector<vector<int>> array = {{}};
    CSCMatrix<int> m1 = from_vector<int>(array);
    CHECK_THROWS_WITH_AS(find_min_valueCSC<int>(m1), "The matrix is empty.", std::exception);
}

TEST_CASE("testing CSC scalar multiplication")
{
    vector<vector<int>> array = {{1, 0, 0}, {4, 5, 6}, {0, 8, 9}};
    CSCMatrix<int> m1 = from_vector<int>(array);
    CSCMatrix<int> m2 = scalar_multiply_matrixCSC<int>(m1, 2);
    vector<vector<int>> scalarMultiplyResultExpected = {{2, 0, 0}, {8, 10, 12}, {0, 16, 18}};
    // check the scalar multiply
    CHECKCSC(m2, scalarMultiplyResultExpected);
}

TEST_CASE("testing CSC scalar multiplication exceptions")
{
    vector<vector<int>> array = {{}};
    CSCMatrix<int> m1 = from_vector<int>(array);
    CHECK_THROWS_WITH_AS(scalar_multiply_matrixCSC<int>(m1, 2), "The matrix is empty.", std::exception);
}

// test for Gaussian Elimination
TEST_CASE("testing CSC Gaussian Elimination 1")
{
    std::vector<std::vector<double>> A1 = {{2, 1, -1},
                                           {-3, -1, 2},
                                           {-2, 1, 2}};
    std::vector<double> b1 = {8, -11, -3};
    bool result1 = gaussian_elimination(A1, b1);
    REQUIRE(result1);
    REQUIRE(Approx(b1[0]).epsilon(0.01) == 2);
    REQUIRE(Approx(b1[1]).epsilon(0.01) == 3);
    REQUIRE(Approx(b1[2]).epsilon(0.01) == -1);
}

TEST_CASE("testing CSC Gaussian Elimination 2")
{
    std::vector<std::vector<double>> A2 = {{1, 2, -1},
                                           {2, 1, -2},
                                           {-3, 1, 1}};
    std::vector<double> b2 = {3, 3, -6};
    bool result2 = gaussian_elimination(A2, b2);
    REQUIRE(result2);
    REQUIRE(Approx(b2[0]).epsilon(0.01) == 1);
    REQUIRE(Approx(b2[1]).epsilon(0.01) == -2);
    REQUIRE(Approx(b2[2]).epsilon(0.01) == -2);
}

TEST_CASE("testing CSC Gaussian Elimination 3")
{
    std::vector<std::vector<double>> A3 = {{1, 2},
                                           {2, 4}};
    std::vector<double> b3 = {3, 6};
    bool result3 = gaussian_elimination(A3, b3);
    REQUIRE(!result3);
}

// test for QR decomposition
TEST("QR Factorization Test 1")
{
    // Create a matrix
    std::vector<std::vector<double>> A = {{1.0, 2.0, 3.0},
                                          {4.0, 5.0, 6.0},
                                          {7.0, 8.0, 7.0},
                                          {4.0, 2.0, 1.0}};

    // Call the QR factorization function
    auto qr = QR_factorization(A);

    // Check that the dimensions of the Q and R matrices are correct
    EXPECT_EQ(qr.Q.size(), A.size());
    EXPECT_EQ(qr.Q[0].size(), A[0].size());
    EXPECT_EQ(qr.R.size(), A[0].size());
    EXPECT_EQ(qr.R[0].size(), A[0].size());

    // Check that Q is orthogonal
    std::vector<std::vector<double>> Q_transpose = transpose(qr.Q);
    std::vector<std::vector<double>> identity = matrix_multiplication(Q_transpose, qr.Q);
    for (size_t i = 0; i < identity.size(); i++)
    {
        for (size_t j = 0; j < identity[0].size(); j++)
        {
            if (i == j)
            {
                EXPECT_NEAR(identity[i][j], 1.0, 1e-6);
            }
            else
            {
                EXPECT_NEAR(identity[i][j], 0.0, 1e-6);
            }
        }
    }

    // Check that Q*R = A
    std::vector<std::vector<double>> Q_times_R = matrix_multiplication(qr.Q, qr.R);
    for (size_t i = 0; i < A.size(); i++)
    {
        for (size_t j = 0; j < A[0].size(); j++)
        {
            EXPECT_NEAR(Q_times_R[i][j], A[i][j], 1e-6);
        }
    }
}

// test for LU factorization
TEST_CASE("LU Factorization test 1")
{
    // Test a 2x2 matrix
    std::vector<std::vector<double>> A1 = {{4, 3}, {6, 3}};
    std::vector<std::vector<double>> expected_L1 = {{1, 0}, {2.0 / 3, 1}};
    std::vector<std::vector<double>> expected_U1 = {{4, 3}, {0, -1}};
    std::vector<int> expected_p1 = {0, 1};

    std::vector<int> p1 = lu_factorization(A1);
    REQUIRE(p1 == expected_p1);
    REQUIRE(A1 == expected_L1);             // The lower triangular part is overwritten with L
    REQUIRE(A1[0][0] == expected_U1[0][0]); // The diagonal and upper triangular part is U
    REQUIRE(A1[0][1] == expected_U1[0][1]);
    REQUIRE(A1[1][0] == expected_U1[1][0]);
    REQUIRE(A1[1][1] == expected_U1[1][1]);
}

TEST_CASE("LU Factorization test 2")
{
    // Test a 3x3 matrix
    std::vector<std::vector<double>> A2 = {{1, 2, 3}, {2, 5, 3}, {1, 0, 8}};
    std::vector<std::vector<double>> expected_L2 = {{1, 0, 0}, {2, 1, 0}, {1, -4, 1}};
    std::vector<std::vector<double>> expected_U2 = {{1, 2, 3}, {0, 1, -3}, {0, 0, 1}};
    std::vector<int> expected_p2 = {0, 1, 2};

    std::vector<int> p2 = lu_factorization(A2);
    REQUIRE(p2 == expected_p2);
    REQUIRE(A2 == expected_L2);             // The lower triangular part is overwritten with L
    REQUIRE(A2[0][0] == expected_U2[0][0]); // The diagonal and upper triangular part is U
    REQUIRE(A2[0][1] == expected_U2[0][1]);
    REQUIRE(A2[0][2] == expected_U2[0][2]);
    REQUIRE(A2[1][1] == expected_U2[1][1]);
    REQUIRE(A2[1][2] == expected_U2[1][2]);
    REQUIRE(A2[2][2] == expected_U2[2][2]);
}

// test for LDL factorization
TEST_CASE("LDL^T Factorization")
{
    std::vector<std::vector<double>> A = {{4, 12, -16}, {12, 37, -43}, {-16, -43, 98}};
    std::vector<double> d;
    bool success = ldlt_factorization(A, d);

    REQUIRE(success == true);

    // Check that LDL^T = A
    const int n = A.size();
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double sum = 0.0;
            for (int k = 0; k < n; k++)
            {
                sum += A[i][k] * d[k] * A[j][k];
            }
            REQUIRE(std::abs(sum - A[i][j]) < 1e-10);
        }
    }
}

TEST_CASE("Gauss-Seidel")
{
    std::vector<std::vector<double>> A = {{4, 1, -1},
                                          {2, 5, 2},
                                          {1, 2, 4}};
    std::vector<double> b = {4, 7, 14};
    std::vector<double> x = {0, 0, 0};

    int max_iter = 100;
    bool success = gauss_seidel(A, b, x, max_iter);

    REQUIRE(success);
    REQUIRE(std::abs(x[0] - 1) < 1e-6);
    REQUIRE(std::abs(x[1] - 1) < 1e-6);
    REQUIRE(std::abs(x[2] - 3) < 1e-6);
}

TEST_CASE("Jacobi Iteration")
{
    const std::vector<std::vector<double>> A = {{4.0, 1.0, 1.0}, {1.0, 4.0, 1.0}, {1.0, 1.0, 4.0}};
    const std::vector<double> b = {6.0, 6.0, 6.0};
    const double tol = 1e-6;
    const int max_iter = 100;

    const std::vector<double> x_expected = {1.0, 1.0, 1.0};
    const std::vector<double> x = jacobi_iteration(A, b, tol, max_iter);

    for (int i = 0; i < x.size(); i++)
    {
        REQUIRE(std::abs(x[i] - x_expected[i]) < tol);
    }
}