#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../functionsCSRParallel.cc"
#include "fstream"
//Basic Unit tests for CSR add, multiply, and transpose
//Use -d to time the tests

/*
Ensure two CSR matrices are the same
*/
void CHECKCSR(CSRMatrix<double> one, CSRMatrix<double> two){
    CHECK(one.numRows == two.numRows);
    CHECK(one.numColumns == two.numColumns);
    for(size_t i = 0; i < one.numRows;i++){
        for(size_t j = 0; j< one.numColumns;j++){
             CHECK_MESSAGE(get_matrixCSR(one, i, j) == get_matrixCSR(two, i, j),"i = "<< i <<", j = "<<j);
        }
    }
}

// TEST_CASE("CSR ADD parallel Correctness") {

//     CSRMatrix<double> m1 = load_fileCSR<double>("../../data/matrices/1138_bus.mtx");
//     CSRMatrix<double> m2 = load_fileCSR<double>("../../data/matrices/1138_bus.mtx");

//     CSRMatrix<double> m3 = add_matrixCSR<double>(m1, m2);
//     CSRMatrix<double> m4 = parallel::add_matrixCSR<double>(m1, m2);
//     CHECKCSR(m3,m4);
// }


TEST_CASE("CSR ADD parallel TIME") {

    CSRMatrix<double> m1 = load_fileCSR<double>("../../data/matrices/1138_bus.mtx");

    CSRMatrix<double> m2 = load_fileCSR<double>("../../data/matrices/1138_bus.mtx");
    //for(size_t i=0 ; i < 1000;i++){
    CSRMatrix<double> m3 = parallel::add_matrixCSR<double>(m1, m2);
    //}
}

TEST_CASE("CSR ADD serial TIME") {

    CSRMatrix<double> m1 = load_fileCSR<double>("../../data/matrices/1138_bus.mtx");

    CSRMatrix<double> m2 = load_fileCSR<double>("../../data/matrices/1138_bus.mtx");
    //for(size_t i=0 ; i < 1000;i++){
    CSRMatrix<double> m3 = add_matrixCSR<double>(m1, m2);
    //}
}

TEST_CASE("CSR max parallel Correctness") {

    CSRMatrix<double> m1 = load_fileCSR<double>("../../data/matrices/1138_bus.mtx");

    double serial = find_max_CSR<double>(m1);
    double parallel = parallel::find_max_CSR<double>(m1);
    CHECK(serial == parallel);
}

TEST_CASE("CSR max parallel TIME") {

    CSRMatrix<double> m1 = load_fileCSR<double>("../../data/matrices/kmer_V1r.mtx");

    double parallel = parallel::find_max_CSR<double>(m1);
    CHECK(parallel == parallel);
}

TEST_CASE("CSR max serial TIME") {

    CSRMatrix<double> m1 = load_fileCSR<double>("../../data/matrices/kmer_V1r.mtx");

    double serial = find_max_CSR<double>(m1);
    CHECK(serial == serial);
}
