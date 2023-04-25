#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../functionsCOOParallel.cc"
#include "../functionsCOO.cc"


void CHECKCOO(COO::COOMatrix<int> mResult, std::vector<std::vector<int>> mCheck) {
    CHECK(mResult.numRows == mCheck.size());
    for(size_t i = 0; i < mResult.numRows;i++){
        CHECK(mResult.numCols == mCheck[i].size());
        for(size_t j = 0; j< mResult.numCols;j++){
             CHECK_MESSAGE(get_matrixCOO(mResult, i, j) == mCheck[i][j],"i = "<< i <<", j = "<<j);
        }
    }
}

TEST_CASE("COO max parallel TIME") {
    COO::COOMatrix<double> m1 = COO::load_fileCOO<double>("../../data/matrices/kmer_V1r.mtx");

    double parallel = COOPar
}

TEST_CASE("COO max serial TIME") {
    COO::COOMatrix<double> m1 = COO::load_fileCOO<double>("../../data/matrices/kmer_V1r.mtx");

    double serial = COO::find_max_COO<double>(m1);
    CHECK(serial == serial);
}

