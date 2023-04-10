#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <tbb/tbb.h>
#include "functionsCSR.cc"

#include <chrono>

namespace parallel {
using namespace std;
/// @brief Adds two compressed spares row(CSR) matrixes together
/// @exception The two matrixes must have the same dimensions
/// @tparam T The type of both matrixes
/// @param m1 The first matrix too add
/// @param m2 The second matrix too add
/// @return m1+m2
template<typename T>
CSRMatrix<T> add_matrixCSR(CSRMatrix<T> m1, CSRMatrix<T> m2){
    if(m1.numRows!= m2.numRows){
        throw std::invalid_argument("The number of rows in the first matrix must match the number of rows in the second matrix.");
    }
    if(m1.numColumns!= m2.numColumns){
        throw std::invalid_argument("The number of columns in the first matrix must match the number of columns in the second matrix.");
    }
    CSRMatrix<T> returnMatrix = tbb::parallel_reduce(tbb::blocked_range<int>(0, m1.numRows), CSRMatrix<T>(),
        [m1,m2](const tbb::blocked_range<int>& r, CSRMatrix<T> v) -> CSRMatrix<T> {
            for (auto i = r.begin(); i < r.end(); i++) {
                size_t a1 = m1.row_ptr.at(i);
                size_t b1 = m1.row_ptr.at(i+1);
                size_t a2 = m2.row_ptr.at(i);
                size_t b2 = m2.row_ptr.at(i+1);
                while(a1 < b1 && a2 < b2){
                    if (m1.col_ind.at(a1) < m2.col_ind.at(a2)) {
                        v.val.push_back(m1.val.at(a1));
                        v.col_ind.push_back(m1.col_ind.at(a1));
                        a1++;
                    } else if (m1.col_ind.at(a1) > m2.col_ind.at(a2)) {
                        v.val.push_back(m2.val.at(a2));
                        v.col_ind.push_back(m2.col_ind.at(a2));
                        a2++;
                    }
                    else if (m1.col_ind.at(a1) == m2.col_ind.at(a2)) {
                        T value = m1.val.at(a1) + m2.val.at(a2);
                        if (value != 0) {
                            v.val.push_back(value);
                            v.col_ind.push_back(m1.col_ind.at(a1));
                        }
                        a1++;
                        a2++;
                    }
                }
                while (a1 < b1) {
                    v.val.push_back(m1.val.at(a1));
                    v.col_ind.push_back(m1.col_ind.at(a1));
                    a1++;
                }
                while (a2 < b2) {
                    v.val.push_back(m2.val.at(a2));
                    v.col_ind.push_back(m2.col_ind.at(a2));
                    a2++;
                }
                v.row_ptr.push_back(v.val.size());
            }
            return v;
        },
        [m1,m2](CSRMatrix<T> v1, CSRMatrix<T> v2) -> CSRMatrix<T> {
            v1.row_ptr.insert(v1.row_ptr.end(),v2.row_ptr.cbegin(),v2.row_ptr.cend());
            v1.col_ind.insert(v1.col_ind.end(),v2.col_ind.cbegin(),v2.col_ind.cend());
            v1.val.insert(v1.val.end(),v2.val.cbegin(),v2.val.cend());
            return v1;
        }
    );
    returnMatrix.numRows = m1.numRows;
    returnMatrix.numColumns = m1.numColumns;
    returnMatrix.row_ptr.insert(returnMatrix.row_ptr.begin(),1,0);
    return returnMatrix;
}

/// @brief Find the max value in a compressed sparse row(CSR) matrix
/// @tparam T The type of the matrix
/// @param m The CSR matrix to find the max value of
/// @return The max value in the matrix
template <typename T>
T find_max_CSR(CSRMatrix<T> m1)
{
    return tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, m1.val.size()),
        m1.val[0],
        [&](const tbb::blocked_range<size_t>& r, T max_value) {
            for (size_t i = r.begin(); i != r.end(); ++i)
            {
                if (m1.val[i] > max_value)
                {
                    max_value = m1.val[i];
                }
            }
            return max_value;
        },
        [](T x, T y) { return std::max(x, y); }
    );
}
}