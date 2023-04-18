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

    
}
