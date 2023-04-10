#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../functionsCSRParallel.cc"
#include "../functionsCSR.cc"
#include "fstream"
//Basic Unit tests for CSR add, multiply, and transpose
//Use -d to time the tests

/*
Takes a CSRMatrix and a dense matrix and makes sure they have the same elements
*/
/// @brief Ensure compressed sparse row matrix is the same as expected dense matrix
/// @param mResult The CSR matrix to compare to
/// @param mCheck The dense matrix to compare to
// void CHECKCSR(CSRMatrix<int> mResult, vector<vector<int> > mCheck){
//     CHECK(mResult.numRows == mCheck.size());
//     for(size_t i = 0; i < mResult.numRows;i++){
//         CHECK(mResult.numColumns == mCheck[i].size());
//         for(size_t j = 0; j< mResult.numColumns;j++){
//              CHECK_MESSAGE(get_matrixCSR(mResult, i, j) == mCheck[i][j],"i = "<< i <<", j = "<<j);
//         }
//     }
// }


TEST_CASE("testing CSR on real matrix") {

    parallel::CSRMatrix<double> m1 = parallel::load_fileCSR<double>("../../data/matrices/1138_bus.mtx");

    parallel::CSRMatrix<double> m2 = parallel::load_fileCSR<double>("../../data/matrices/1138_bus.mtx");
    //for(size_t i=0 ; i < 1000;i++){
    parallel::CSRMatrix<double> m3 = parallel::add_matrixCSR<double>(m1, m2);
    //}


}
