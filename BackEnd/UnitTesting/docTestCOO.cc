#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../functionsCOO.cc"

#include <vector>

using namespace std;

void CHECKCOO(COO::COOMatrix<int> mResult, std::vector<std::vector<int>> mCheck) {
    CHECK(mResult.numRows == mCheck.size());
    for(size_t i = 0; i < mResult.numRows;i++){
        CHECK(mResult.numCols == mCheck[i].size());
        for(size_t j = 0; j< mResult.numCols;j++){
             CHECK_MESSAGE(get_matrixCOO(mResult, i, j) == mCheck[i][j],"i = "<< i <<", j = "<<j);
        }
    }
}

TEST_CASE("Testing COO Add") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};

    vector<vector<int> > array2 = {{0,2,3},{0,-5,0},{7,8,0}};

    COO::COOMatrix<int> m1 = COO::from_vector<int>(array);
    CHECKCOO(m1, array);

    COO::COOMatrix<int> m2 = COO::from_vector<int>(array2);
    CHECKCOO(m2, array2);

    COO::COOMatrix<int> m3 = COO::add_matrixCOO(m1, m2);
    vector<vector<int>> addResultExpected = {{1,2,3}, {4,0,6}, {7,16,9}};

    CHECKCOO(m3,addResultExpected);
}

TEST_CASE("COO Add Exceptions") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};
    vector<vector<int> > array2 = {{0,2,3},{0,-5,0}};

    COO::COOMatrix<int> m1 = COO::from_vector<int>(array);
    COO::COOMatrix<int> m2 = COO::from_vector<int>(array2);

    CHECK_THROWS_WITH_AS(COO::add_matrixCOO<int>(m1, m2),"The number of rows in the first matrix must match the number of rows in the second matrix.",std::exception);

     vector<vector<int> > array3 = {{1,0},{4,5},{0,8}};
     COO::COOMatrix<int> m3 = COO::from_vector<int>(array3);

    CHECK_THROWS_WITH_AS(COO::add_matrixCOO<int>(m1, m3),"The number of columns in the first matrix must match the number of columns in the second matrix.",std::exception);
}

TEST_CASE("testing COO transpose") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};


    COO::COOMatrix<int> m1 = COO::from_vector(array);

    COO::COOMatrix<int> m2 = COO::transpose_matrixCOO(m1);
    vector<vector<int> > transposeResultExpected = {{1,4,0},{0,5,8},{0,6,9}};
    //check the transpose
    CHECKCOO(m2,transposeResultExpected);
}

// TEST_CASE("testing COO scalar add") {

// }

TEST_CASE("testing COO multiply") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}}; 
    vector<vector<int> > array2 = {{0,2,3},{0,-5,0},{7,8,0}};

    COO::COOMatrix<int> m1 = COO::from_vector<int>(array);
    COO::COOMatrix<int> m2 = COO::from_vector<int>(array2);
    COO::COOMatrix<int> m3 = COO::multiply_matrixCOO<int>(m1, m2);
    vector<vector<int> > multiplyResultExpected = {{0,2,3},{42,31,12},{63,32,0}};
    //check the multiply
    CHECKCOO(m3,multiplyResultExpected);
}

TEST_CASE("COO multiply zero matrix") {
    vector<vector<int> > array = {{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    vector<vector<int> > array2 = {{0,1},{0,-5},{7,8},{56,76}};

    COO::COOMatrix<int> m1 = COO::from_vector<int>(array);
    COO::COOMatrix<int> m2 = COO::from_vector<int>(array2);
    COO::COOMatrix<int> m3 = COO::multiply_matrixCOO<int>(m1, m2);
    vector<vector<int> > multiplyResultExpected = {{0,0},{0,0},{0,0}};
    //check the multiply
    CHECKCOO(m3,multiplyResultExpected);
}

TEST_CASE("COO multiply zero row first matrix") {
    vector<vector<int> > array = {{0,5,0,3},{0,0,0,0},{0,7,90,0}};
    vector<vector<int> > array2 = {{0,1},{0,-5},{7,8},{56,76}};

    COO::COOMatrix<int> m1 = COO::from_vector<int>(array);
    COO::COOMatrix<int> m2 = COO::from_vector<int>(array2);
    COO::COOMatrix<int> m3 = COO::multiply_matrixCOO<int>(m1, m2);
    vector<vector<int> > multiplyResultExpected = {{168,203},{0,0},{630,685}};
    //check the multiply
    CHECKCOO(m3,multiplyResultExpected);
}

TEST_CASE("COO multiply zero row second matrix") {
    vector<vector<int> > array = {{0,5,0,3},{0,6,7,0},{0,7,90,0}};
    vector<vector<int> > array2 = {{0,0},{0,-5},{0,0},{56,76}};

    COO::COOMatrix<int> m1 = COO::from_vector<int>(array);
    COO::COOMatrix<int> m2 = COO::from_vector<int>(array2);
    COO::COOMatrix<int> m3 = COO::multiply_matrixCOO<int>(m1, m2);
    vector<vector<int> > multiplyResultExpected = {{168,203},{0,-30},{0,-35}};
    //check the multiply
    CHECKCOO(m3,multiplyResultExpected);
}

TEST_CASE("COO multiply zero row both matrix") {
    vector<vector<int> > array = {{0,5,0,3},{0,0,0,0},{0,7,90,0}};
    vector<vector<int> > array2 = {{0,0},{0,-5},{0,0},{56,76}};

    COO::COOMatrix<int> m1 = COO::from_vector<int>(array);
    COO::COOMatrix<int> m2 = COO::from_vector<int>(array2);
    COO::COOMatrix<int> m3 = COO::multiply_matrixCOO<int>(m1, m2);
    vector<vector<int> > multiplyResultExpected = {{168,203},{0,0},{0,-35}};
    //check the multiply
    CHECKCOO(m3,multiplyResultExpected);
}

TEST_CASE("COO multiply Exceptions") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};
    vector<vector<int> > array2 = {{0,2},{0,-5}};

    COO::COOMatrix<int> m1 = COO::from_vector<int>(array);
    COO::COOMatrix<int> m2 = COO::from_vector<int>(array2);

    CHECK_THROWS_WITH_AS(COO::multiply_matrixCOO<int>(m1, m2),"The number of columns in the first matrix must match the number of rows in the second matrix.",std::exception);
}


