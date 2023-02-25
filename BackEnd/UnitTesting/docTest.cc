#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../functionsCSR.cc"
//This file just defines are Doctest config


TEST_CASE("testing CSR") {
    vector<vector<int> > array = vector<vector<int> >(3, vector<int>(3));
    array[0][0] = 1;
    array[1][0] = 4;
    array[1][1] = 5;
    array[1][2] = 6;
    array[2][1] = 8;
    array[2][2] = 9;

    vector<vector<int> > array2 = vector<vector<int> >(3, vector<int>(3));
    array2[0][1] = 2;
    array2[0][2] = 3;
    array2[1][1] = -5;
    array2[2][0] = 7;
    array2[2][1] = 8;

    CSRMatrix<int> m1 = from_vector<int>(array);
    //check m1
    CHECK(get_matrixCSR(m1, 0, 0) == 1);
    CHECK(get_matrixCSR(m1, 1, 0) == 4);
    CHECK(get_matrixCSR(m1, 2, 1) == 8);

    CSRMatrix<int> m2 = from_vector<int>(array2);
    //check m2
    CHECK(get_matrixCSR(m2, 0, 1) == 2);
    CHECK(get_matrixCSR(m2, 1, 1) == -5);
    CHECK(get_matrixCSR(m2, 2, 1) == 8);
    CSRMatrix<int> m3 = add_matrixCSR<int>(m1, m2);
    //check the addition
    CHECK(get_matrixCSR(m3, 0, 0) == 1);
    CHECK(get_matrixCSR(m3, 0, 1) == 2);
    CHECK(get_matrixCSR(m3, 0, 2) == 3);
    CHECK(get_matrixCSR(m3, 1, 0) == 4);
    CHECK(get_matrixCSR(m3, 1, 1) == 0);
    CHECK(get_matrixCSR(m3, 1, 2) == 6);
    CHECK(get_matrixCSR(m3, 2, 0) == 7);
    CHECK(get_matrixCSR(m3, 2, 1) == 16);
    CHECK(get_matrixCSR(m3, 2, 2) == 9);
}