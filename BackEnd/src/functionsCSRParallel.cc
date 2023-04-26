#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <tbb/tbb.h>
#include "functionsCSR.cc"
#include "functions.cc"

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
        tbb::blocked_range<int>(0, m1.val.size()),
        m1.val[0],
        [&](const tbb::blocked_range<int>& r, T max_value) {
            cerr << "hi" << m1.val.size() << endl;
            for (int i = r.begin(); i != r.end(); ++i)
            {
                if (m1.val[i] > max_value)
                {
                    max_value = m1.val[i];
                }
            }
            return max_value;
        },
        [](T x, T y) { return std::max(x, y); }
    );}
//     tbb::parallel_for( tbb::blocked_range<int>(0, 100000000), [&](tbb::blocked_range<int> r){
//         cerr << "hi" << std::endl;
// 			for (int i = r.begin(); i != r.end(); ++i)
//             {
//                m1.val[i%m1.val.size()] = 5 *i;
//             }
// });
//     return true;
// }
// gaussian elimination with partial pivoting
// returns true if successful, false if A is singular
// Modifies both A and b to store the results
bool gaussian_elimination(std::vector<std::vector<double> > &A) {
    // Iterate over each row in the matrix
    float pivot;
    for(size_t i = 0; i < A.size() - 1; i++){
        // Pivot will be the diagonal
        pivot = A[i][i];

        // Iterate of the remaining row elements
        for(size_t j = i + 1; j < A[0].size(); j++){
            // Divide by the pivot
            A[i][j] /= pivot;
        }

        // Do direct assignment for trivial case (self-divide)
        A[i][i] = 1.0;

        // Eliminate ith element from the jth row
        float scale;
        for(size_t j = i + 1; j < A.size(); j++){
            // Factor we will use to scale subtraction by
            scale = A[j][i];

            // Iterate over the remaining columns
            for(size_t k = i + 1; k < A.size(); k++){
                A[j][k] -= A[i][k] * scale;
            }

            // Do direct assignment for trivial case (eliminate position)
            A[j][i] = 0;
        }
    }
    A[A.size()-1][A[0].size()-1] = 1;

    return true;
}

// gaussian elimination with partial pivoting
// returns true if successful, false if A is singular
// Modifies both A and b to store the results
bool gaussian_elimination_parallel(std::vector<std::vector<double> > &A) {
     // Iterate over each row in the matrix
    float pivot;
    for(size_t i = 0; i < A.size() - 1; i++){
        // Pivot will be the diagonal
        pivot = A[i][i];

        // Iterate of the remaining row elements
        tbb::parallel_for( tbb::blocked_range<size_t>(i+1, A[0].size()), [&](tbb::blocked_range<size_t> r){
            for(size_t j = r.begin(); j < r.end(); j++){
                A[i][j] /= pivot;
            }
        });
        // for(size_t j = i + 1; j < A[0].size(); j++){
        //     // Divide by the pivot
        //     A[i][j] /= pivot;
        // }

        // Do direct assignment for trivial case (self-divide)
        A[i][i] = 1.0;

        // Eliminate ith element from the jth row
        tbb::parallel_for( tbb::blocked_range<size_t>(i+1, A.size()), [&](tbb::blocked_range<size_t> r){
            float scale;
            for(size_t j = r.begin(); j < r.end(); j++){
                // Factor we will use to scale subtraction by
                scale = A[j][i];

                // Iterate over the remaining columns
                for(size_t k = i + 1; k < A.size(); k++){
                    A[j][k] -= A[i][k] * scale;
                }

                // Do direct assignment for trivial case (eliminate position)
                A[j][i] = 0;
            }
        });
        // float scale;
        // for(size_t j = i + 1; j < A.size(); j++){
        //     // Factor we will use to scale subtraction by
        //     scale = A[j][i];

        //     // Iterate over the remaining columns
        //     for(size_t k = i + 1; k < A.size(); k++){
        //         A[j][k] -= A[i][k] * scale;
        //     }

        //     // Do direct assignment for trivial case (eliminate position)
        //     A[j][i] = 0;
        // }
    }
    A[A.size()-1][A[0].size()-1] = 1;

    return true;
}
// gaussian elimination with partial pivoting
// returns true if successful, false if A is singular
// Modifies both A and b to store the results
bool gaussian_elimination_swap(std::vector<std::vector<double> > &A) {
    const size_t n = A.size();

    for (size_t i = 0; i < n; i++) {

        // find pivot row and swap
        size_t max = i;
        for (size_t k = i + 1; k < n; k++)
            if (abs(A[k][i]) > abs(A[max][i]))
                max = k;
        std::swap(A[i], A[max]);

        // singular or nearly singular
        if (abs(A[i][i]) <= 1e-10)
            return false;

        // pivot within A and b
        for (size_t k = i + 1; k < n; k++) {
            double t = A[k][i] / A[i][i];
            for (size_t j = i; j < n; j++) {
                A[k][j] -= A[i][j] * t;
            }
        }
    }
    return true;
}

}
#include <chrono>
class timer {
public:
    std::chrono::time_point<std::chrono::high_resolution_clock> lastTime;
    timer() : lastTime(std::chrono::high_resolution_clock::now()) {}
    inline double elapsed() {
        std::chrono::time_point<std::chrono::high_resolution_clock> thisTime=std::chrono::high_resolution_clock::now();
        double deltaTime = std::chrono::duration<double>(thisTime-lastTime).count();
        lastTime = thisTime;
        return deltaTime;
    }
};
#include <iomanip>
void printMatrix(vector<vector<double>> &A) {
    // Set the width of each output element to 8 characters
    const int width = 9;
    // Loop over each row of the matrix
    for (size_t i = 0; i < A.size(); i++) {
        // Loop over each element in the row
        for (size_t j = 0; j < A[i].size(); j++) {
            // Set the output width for the element
            cout << setw(width) << A[i][j] << " ";
        }
        // End the row with a newline
        cout << endl;
    }
}

// int main() {
//     timer stopwatch;
//     std::vector<double> parallel;
//     std::vector<double> serial;
//     for(size_t i = 1 ; i <16;i++){
//         std::vector<std::vector<double> > A = generate_random_matrix(1500,1500,1,10000);
//         std::vector<std::vector<double> > B = generate_random_matrix(1500,1500,1,10000);
//         tbb::task_arena arena(i);
// 	    	arena.execute([&]() {
//             stopwatch.elapsed();
//             parallel::gaussian_elimination_parallel(A);
//             parallel.push_back(stopwatch.elapsed());
//         });
        
        

//         stopwatch.elapsed();
//         parallel::gaussian_elimination(B);
//         serial.push_back(stopwatch.elapsed());
//     }
//     for(auto val : parallel){
//         cerr<< val << ",";
//     }
//     cerr<< endl;
//     for(auto val : serial){
//         cerr<< val << ",";
//     }
//     cerr<< endl;
    


//     return 0;
// }

int main() {
    timer stopwatch;
    std::vector<double> parallel;
    std::vector<double> serial;
    std::vector<int> size = {500,1000,1500,2000,2500,3000};
    for(size_t i = 0 ; i <size.size();i++){
        std::vector<std::vector<double> > A = generate_random_matrix(size[i],size[i],1,10000);
        std::vector<std::vector<double> > B = generate_random_matrix(size[i],size[i],1,10000);
        stopwatch.elapsed();
        parallel::gaussian_elimination_swap(A);
        parallel.push_back(stopwatch.elapsed());
        stopwatch.elapsed();
        parallel::gaussian_elimination(B);
        serial.push_back(stopwatch.elapsed());
    }
    for(auto val : parallel){
        cerr<< val << ",";
    }
    cerr<< endl;
    for(auto val : serial){
        cerr<< val << ",";
    }
    cerr<< endl;
    


    return 0;
}