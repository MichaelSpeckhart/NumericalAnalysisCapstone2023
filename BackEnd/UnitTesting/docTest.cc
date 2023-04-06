#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../functionsCSRParallel.cc"
#include "../functions.cc"
#include "fstream"
//Basic Unit tests for CSR add, multiply, and transpose
//Use -d to time the tests


/// Test the CSRMatrix class
/*
Takes a CSRMatrix and a dense matrix and makes sure they have the same elements
*/
/// @brief Ensure compressed sparse row matrix is the same as expected dense matrix
/// @param mResult The CSR matrix to compare to
/// @param mCheck The dense matrix to compare to
void CHECKCSR(CSRMatrix<int> mResult, vector<vector<int> > mCheck){
    CHECK(mResult.numRows == mCheck.size());
    for(size_t i = 0; i < mResult.numRows;i++){
        CHECK(mResult.numColumns == mCheck[i].size());
        for(size_t j = 0; j< mResult.numColumns;j++){
             CHECK_MESSAGE(get_matrixCSR(mResult, i, j) == mCheck[i][j],"i = "<< i <<", j = "<<j);
        }
    }
}

TEST_CASE("testing CSR Add") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};

    vector<vector<int> > array2 = {{0,2,3},{0,-5,0},{7,8,0}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    //check m1
    CHECKCSR(m1,array);

    CSRMatrix<int> m2 = from_vector<int>(array2);
    //check m2
    CHECKCSR(m2,array2);

    CSRMatrix<int> m3 = add_matrixCSR<int>(m1, m2);
    vector<vector<int> > addResultExpected = {{1,2,3},{4,0,6},{7,16,9}};
    //check the addition
    CHECKCSR(m3,addResultExpected);
}

TEST_CASE("CSR Add Exceptions") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};
    vector<vector<int> > array2 = {{0,2,3},{0,-5,0}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);

    CHECK_THROWS_WITH_AS(add_matrixCSR<int>(m1, m2),"The number of rows in the first matrix must match the number of rows in the second matrix.",std::exception);

     vector<vector<int> > array3 = {{1,0},{4,5},{0,8}};
     CSRMatrix<int> m3 = from_vector<int>(array3);

    CHECK_THROWS_WITH_AS(add_matrixCSR<int>(m1, m3),"The number of columns in the first matrix must match the number of columns in the second matrix.",std::exception);
}

TEST_CASE("testing CSR transpose") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};


    CSRMatrix<int> m1 = from_vector<int>(array);

    CSRMatrix<int> m2 = transpose_matrixCSR<int>(m1);
    vector<vector<int> > transposeResultExpected = {{1,4,0},{0,5,8},{0,6,9}};
    //check the transpose
    CHECKCSR(m2,transposeResultExpected);
}

TEST_CASE("testing CSR multiply") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};
    vector<vector<int> > array2 = {{0,2,3},{0,-5,0},{7,8,0}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = multiply_matrixCSR<int>(m1, m2);
    vector<vector<int> > multiplyResultExpected = {{0,2,3},{42,31,12},{63,32,0}};
    //check the multiply
    CHECKCSR(m3,multiplyResultExpected);
}

TEST_CASE("CSR multiply zero matrix") {
    vector<vector<int> > array = {{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    vector<vector<int> > array2 = {{0,1},{0,-5},{7,8},{56,76}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = multiply_matrixCSR<int>(m1, m2);
    vector<vector<int> > multiplyResultExpected = {{0,0},{0,0},{0,0}};
    //check the multiply
    CHECKCSR(m3,multiplyResultExpected);
}

TEST_CASE("CSR multiply zero row first matrix") {
    vector<vector<int> > array = {{0,5,0,3},{0,0,0,0},{0,7,90,0}};
    vector<vector<int> > array2 = {{0,1},{0,-5},{7,8},{56,76}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = multiply_matrixCSR<int>(m1, m2);
    vector<vector<int> > multiplyResultExpected = {{168,203},{0,0},{630,685}};
    //check the multiply
    CHECKCSR(m3,multiplyResultExpected);
}

TEST_CASE("CSR multiply zero row second matrix") {
    vector<vector<int> > array = {{0,5,0,3},{0,6,7,0},{0,7,90,0}};
    vector<vector<int> > array2 = {{0,0},{0,-5},{0,0},{56,76}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = multiply_matrixCSR<int>(m1, m2);
    vector<vector<int> > multiplyResultExpected = {{168,203},{0,-30},{0,-35}};
    //check the multiply
    CHECKCSR(m3,multiplyResultExpected);
}

TEST_CASE("CSR multiply zero row both matrix") {
    vector<vector<int> > array = {{0,5,0,3},{0,0,0,0},{0,7,90,0}};
    vector<vector<int> > array2 = {{0,0},{0,-5},{0,0},{56,76}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = multiply_matrixCSR<int>(m1, m2);
    vector<vector<int> > multiplyResultExpected = {{168,203},{0,0},{0,-35}};
    //check the multiply
    CHECKCSR(m3,multiplyResultExpected);
}

TEST_CASE("CSR multiply Exceptions") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};
    vector<vector<int> > array2 = {{0,2},{0,-5}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);

    CHECK_THROWS_WITH_AS(multiply_matrixCSR<int>(m1, m2),"The number of columns in the first matrix must match the number of rows in the second matrix.",std::exception);
}

TEST_CASE("testing CSR multiply scalar") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = multiply_matrixCSR<int>(m1, 2);
    vector<vector<int> > multiplyResultExpected = {{2,0,0},{8,10,12},{0,16,18}};
    //check the multiply
    CHECKCSR(m2,multiplyResultExpected);
}

TEST_CASE("testing CSR multiply scalar zero matrix") {
    vector<vector<int> > array = {{0,0,0},{0,0,0},{0,0,0}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = multiply_matrixCSR<int>(m1, 2);
    vector<vector<int> > multiplyResultExpected = {{0,0,0},{0,0,0},{0,0,0}};
    //check the multiply
    CHECKCSR(m2,multiplyResultExpected);
}

TEST_CASE("testing CSR multiply scalar zero row") {
    vector<vector<int> > array = {{0,0,0},{0,5,0},{0,0,0}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = multiply_matrixCSR<int>(m1, 2);
    vector<vector<int> > multiplyResultExpected = {{0,0,0},{0,10,0},{0,0,0}};
    //check the multiply
    CHECKCSR(m2,multiplyResultExpected);
}

TEST_CASE("testing CSR multiply scalar zero column") {
    vector<vector<int> > array = {{0,0,0},{0,5,0},{0,0,0}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = multiply_matrixCSR<int>(m1, 2);
    vector<vector<int> > multiplyResultExpected = {{0,0,0},{0,10,0},{0,0,0}};
    //check the multiply
    CHECKCSR(m2,multiplyResultExpected);
}

TEST_CASE("testing CSR multiply scalar zero row and column") {
    vector<vector<int> > array = {{0,0,0},{0,0,0},{0,0,0}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = multiply_matrixCSR<int>(m1, 2);
    vector<vector<int> > multiplyResultExpected = {{0,0,0},{0,0,0},{0,0,0}};
    //check the multiply
    CHECKCSR(m2,multiplyResultExpected);
}

TEST_CASE("testing CSR subtract") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};
    vector<vector<int> > array2 = {{0,2},{0,-5},{7,8}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = subtract_matrixCSR<int>(m1, m2);
    vector<vector<int> > subtractResultExpected = {{1,-2,0},{4,10,6},{-7,0,1}};
    //check the subtract
    CHECKCSR(m3,subtractResultExpected);
}

TEST_CASE("testing CSR subtract zero matrix") {
    vector<vector<int> > array = {{0,0,0},{0,0,0},{0,0,0}};
    vector<vector<int> > array2 = {{0,0,0},{0,0,0},{0,0,0}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = subtract_matrixCSR<int>(m1, m2);
    vector<vector<int> > subtractResultExpected = {{0,0,0},{0,0,0},{0,0,0}};
    //check the subtract
    CHECKCSR(m3,subtractResultExpected);
}

TEST_CASE("CSR Subtract Exceptions") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};
    vector<vector<int> > array2 = {{0,2},{0,-5}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);

    CHECK_THROWS_WITH_AS(subtract_matrixCSR<int>(m1, m2),"The number of columns in the first matrix must match the number of rows in the second matrix.",std::exception);
}

TEST_CASE("testing CSR subtract zero row") {
    vector<vector<int> > array = {{0,0,0},{0,5,0},{0,0,0}};
    vector<vector<int> > array2 = {{0,0,0},{0,0,0},{0,0,0}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = subtract_matrixCSR<int>(m1, m2);
    vector<vector<int> > subtractResultExpected = {{0,0,0},{0,5,0},{0,0,0}};
    //check the subtract
    CHECKCSR(m3,subtractResultExpected);
}

TEST_CASE("CSR Subtract Exceptions") {
    //number of rows in m1 must match number of rows in m2
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};
    vector<vector<int> > array2 = {{0,2},{0,-5}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CHECK_THROWS_WITH_AS(subtract_matrixCSR<int>(m1, m2),"The number of rows in the first matrix must match the number of rows in the second matrix.",std::exception);
}

TEST_CASE("testing CSR subtract zero column") {
    vector<vector<int> > array = {{0,0,0},{0,5,0},{0,0,0}};
    vector<vector<int> > array2 = {{0,0,0},{0,0,0},{0,0,0}};

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = subtract_matrixCSR<int>(m1, m2);
    vector<vector<int> > subtractResultExpected = {{0,0,0},{0,5,0},{0,0,0}};
    //check the subtract
    CHECKCSR(m3,subtractResultExpected);
}

TEST_CASE("Find min value") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    int min = find_min_value<int>(m1);
    CHECK(min == 0);
}

TEST_CASE("Find min value zero matrix") {
    vector<vector<int> > array = {{0,0,0},{0,0,0},{0,0,0}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    int min = find_min_value<int>(m1);
    CHECK(min == 0);
}

TEST_CASE("Find min value zero row") {
    vector<vector<int> > array = {{0,0,0},{0,5,0},{0,0,0}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    int min = find_min_value<int>(m1);
    CHECK(min == 0);
}

TEST_CASE("Find min value zero column") {
    vector<vector<int> > array = {{0,0,0},{0,5,0},{0,0,0}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    int min = find_min_value<int>(m1);
    CHECK(min == 0);
}

TEST_CASE("Find max value") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    int max = find_max_value<int>(m1);
    CHECK(max == 9);
}

TEST_CASE("Find max value zero matrix") {
    vector<vector<int> > array = {{0,0,0},{0,0,0},{0,0,0}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    int max = find_max_value<int>(m1);
    CHECK(max == 0);
}

TEST_CASE("Find max value zero row") {
    vector<vector<int> > array = {{0,0,0},{0,5,0},{0,0,0}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    int max = find_max_value<int>(m1);
    CHECK(max == 5);
}

TEST_CASE("Find max value zero column") {
    vector<vector<int> > array = {{0,0,0},{0,5,0},{0,0,0}};
    CSRMatrix<int> m1 = from_vector<int>(array);
    int max = find_max_value<int>(m1);
    CHECK(max == 5);
}


/// Test the CSCMatrix class


