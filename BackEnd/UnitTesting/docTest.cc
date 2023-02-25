#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../functionsCSR.cc"
//This file just defines are Doctest config


void CHECKCSR(CSRMatrix<int> mResult, vector<vector<int> > mCheck){
    CHECK(mResult.numRows == mCheck.size());
    for(size_t i = 0; i < mResult.numRows;i++){
        CHECK(mResult.numColumns == mCheck[i].size());
        for(size_t j = 0; j< mResult.numColumns;j++){
             CHECK(get_matrixCSR(mResult, i, j) == mCheck[i][j]);
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

TEST_CASE("testing CSR transpose") {
    vector<vector<int> > array = {{1,0,0},{4,5,6},{0,8,9}};


    CSRMatrix<int> m1 = from_vector<int>(array);

    CSRMatrix<int> m2 = transpose_matrixCSR<int>(m1);
    vector<vector<int> > transposeResultExpected = {{1,4,0},{0,5,8},{0,6,9}};
    //check the transpose
    CHECKCSR(m2,transposeResultExpected);
}