#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "../functionsCSRParallel.cc"
#include "../functionsParallel.cc"
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
void CHECKCSRFast(CSRMatrix<double> one, CSRMatrix<double> two){
    CHECK(one.numRows == two.numRows);
    CHECK(one.numColumns == two.numColumns);
    CHECK(one.row_ptr.size() == two.row_ptr.size());
    CHECK(one.col_ind.size() == two.col_ind.size());
    CHECK(one.val.size() == two.val.size());
    for(size_t i = 0 ; i <one.row_ptr.size();i++){
        CHECK_MESSAGE(one.row_ptr[i] == two.row_ptr[i],"i = "<< i);
    }
    for(size_t i = 0 ; i <one.col_ind.size();i++){
         CHECK_MESSAGE(one.col_ind[i] == two.col_ind[i],"i = "<< i);
    }
    for(size_t i = 0 ; i <one.val.size();i++){
        CHECK_MESSAGE(one.val[i] == two.val[i],"i = "<< i);
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

// TEST_CASE("CSR ADD parallel Correctness") {

//     CSRMatrix<double> m1 = load_fileCSR<double>("../../../data/matrices/1138_bus.mtx");
//     CSRMatrix<double> m2 = load_fileCSR<double>("../../../data/matrices/1138_bus.mtx");

//     CSRMatrix<double> m3 = add_matrixCSR<double>(m1, m2);
//     CSRMatrix<double> m4 = parallel::add_matrixCSR<double>(m1, m2);
//     CHECKCSR(m3,m4);
// }


// TEST_CASE("CSR ADD parallel TIME") {

//     CSRMatrix<double> m1 = load_fileCSR<double>("../../../data/matrices/1138_bus.mtx");

//     CSRMatrix<double> m2 = load_fileCSR<double>("../../../data/matrices/1138_bus.mtx");
//     //for(size_t i=0 ; i < 1000;i++){
//     CSRMatrix<double> m3 = parallel::add_matrixCSR<double>(m1, m2);
//     //}
// }

// TEST_CASE("CSR ADD serial TIME") {

//     CSRMatrix<double> m1 = load_fileCSR<double>("../../../data/matrices/1138_bus.mtx");

//     CSRMatrix<double> m2 = load_fileCSR<double>("../../../data/matrices/1138_bus.mtx");
//     //for(size_t i=0 ; i < 1000;i++){
//     CSRMatrix<double> m3 = add_matrixCSR<double>(m1, m2);
//     //}
// }

// TEST_CASE("CSR max parallel Correctness") {

//     CSRMatrix<double> m1 = load_fileCSR<double>("../../../data/matrices/1138_bus.mtx");

//     double serial = find_max_CSR<double>(m1);
//     double parallel = parallel::find_max_CSR<double>(m1);
//     CHECK(serial == parallel);
// }

// TEST_CASE("CSR max parallel TIME") {

//     CSRMatrix<double> m1 = load_fileCSR<double>("../../../data/matrices/kmer_V1r.mtx");

//     double parallel = parallel::find_max_CSR<double>(m1);
//     CHECK(parallel ==parallel);
// }

// TEST_CASE("CSR max serial TIME") {

//     CSRMatrix<double> m1 = load_fileCSR<double>("../../../data/matrices/kmer_V1r.mtx");
//     double serial = find_max_CSR<double>(m1);
//     CHECK(serial == serial);
// }

// TEST_CASE("CSR multiply parallel Correctness") {
//     CSRMatrix<double> m1 = load_fileCSR<double>("../../../data/matrices/TSOPF_RS_b39_c30.mtx");
//     //CSRMatrix<double> m1 = load_fileCSR<double>("../../../data/matrices/1138_bus.mtx");
//     //CSRMatrix<double> m2 = transpose_matrixCSR<double>(m1);

//     CSRMatrix<double> m3 = multiply_matrixCSR<double>(m1, m1);
//     CSRMatrix<double> m4 = parallel::multiply_matrixCSR<double>(m1, m1);
//     CHECKCSRFast(m3,m4);
//     // CHECK(m3.numRows == m4.numRows);
//     // CHECK(m3.numColumns == m4.numColumns);
//     // for(size_t i = 0 ; i <m3.row_ptr.size();i++){
//     //     CHECK(m3.row_ptr[i] == m4.row_ptr[i]);
//     // }
//     // for(size_t i = 0 ; i <m3.col_ind.size();i++){
//     //     CHECK(m3.col_ind[i] == m4.col_ind[i]);
//     // }
//     // for(size_t i = 0 ; i <m3.val.size();i++){
//     //     CHECK(m3.val[i] == m4.val[i]);
//     // }
// }


TEST_CASE("Jacobi Iteration CSR 0")
{
    std::vector<std::vector<double>> A = {{4.0, 1.0, 1.0}, {1.0, 4.0, 1.0}, {1.0, 1.0, 4.0}};
    auto CSR_A = from_vector_CSR<double>(A);
    std::vector<double> b = {6.0, 6.0, 6.0};
    const double tol = 1e-6;
    const int max_iter = 100;

    std::vector<double> x_expected = {1.0, 1.0, 1.0};
    std::vector<double> x = parallel::jacobi_method_CSR<double>(CSR_A, b,tol, max_iter);

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
    std::vector<double> x = parallel::jacobi_method_CSR(CSR_A, b, tol, max_iter);

    CHECK_VECTOR_EQ(x, x_expected, tol);
}

TEST_CASE("Gauss Seidel Iteration CSR 0")
{
    std::vector<std::vector<double>> A = {{4.0, 1.0, 1.0}, {1.0, 4.0, 1.0}, {1.0, 1.0, 4.0}};
    auto CSR_A = from_vector_CSR<double>(A);
    std::vector<double> b = {6.0, 6.0, 6.0};
    const double tol = 1e-6;
    const int max_iter = 100;

    std::vector<double> x_expected = {1.0, 1.0, 1.0};
    std::vector<double> x = parallel::gauss_sidel_CSR<double>(CSR_A, b,tol, max_iter);

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
    std::vector<double> x = parallel::gauss_sidel_CSR(CSR_A, b, tol, max_iter);

    CHECK_VECTOR_EQ(x, x_expected, tol);
}