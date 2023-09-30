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

namespace parallel {
using namespace std;

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

template<typename T>
//to store the column and value vectors in the CSR format
class VectorPair{
    public:
    vector<size_t> vec1;
    vector<T> vec2;

    VectorPair() : vec1(), vec2() {}
};

/// @brief Adds two compressed spares row(CSR) matrixes together
/// @exception The two matrixes must have the same dimensions
/// @tparam T The type of both matrixes
/// @param m1 The first matrix too add
/// @param m2 The second matrix too add
/// @return m1+m2
template <typename T>
CSRMatrix<T> add_matrixCSR(CSRMatrix<T> m1, CSRMatrix<T> m2)
{
    if (m1.numRows != m2.numRows)
    {
        throw std::invalid_argument("The number of rows in the first matrix must match the number of rows in the second matrix.");
    }
    if (m1.numColumns != m2.numColumns)
    {
        throw std::invalid_argument("The number of columns in the first matrix must match the number of columns in the second matrix.");
    }
    CSRMatrix<T> returnMatrix;
    returnMatrix.numRows = m1.numRows;
    returnMatrix.numColumns = m1.numColumns;
    returnMatrix.row_ptr.push_back(0);
    //set the processor count to the number of available number of threads
    const auto processor_count = std::thread::hardware_concurrency();
    //vector of vector pairs to store the results of each thread
    std::vector<vector<VectorPair<T>>> result;
    //resize to ensure each thread has a place to put it's results
    result.resize(processor_count);
    //the vector of threads
  	std::vector<std::thread> threadPool;
	threadPool.reserve(processor_count);
    //lambda to do the add
    auto do_work = [&](int k){
        vector<VectorPair<T>> v;
        //increments by processor count and starts at k to divide work evenly
		for (size_t i = k; i < m1.numRows; i += processor_count){
        VectorPair<T> row;
        size_t a1 = m1.row_ptr.at(i);
        size_t b1 = m1.row_ptr.at(i + 1);
        size_t a2 = m2.row_ptr.at(i);
        size_t b2 = m2.row_ptr.at(i + 1);
        while (a1 < b1 && a2 < b2)
        {
            if (m1.col_ind.at(a1) < m2.col_ind.at(a2))
            {
                row.vec1.push_back(m1.col_ind.at(a1));
                row.vec2.push_back(m1.val.at(a1));
                a1++;
            }
            else if (m1.col_ind.at(a1) > m2.col_ind.at(a2))
            {
                row.vec1.push_back(m2.col_ind.at(a2));
                row.vec2.push_back(m2.val.at(a2));
                a2++;
            }
            else if (m1.col_ind.at(a1) == m2.col_ind.at(a2))
            {
                T value = m1.val.at(a1) + m2.val.at(a2);
                if (value != 0)
                {
                    row.vec1.push_back(m1.col_ind.at(a1));
                    row.vec2.push_back(value);
                }
                a1++;
                a2++;
            }
        }
        while (a1 < b1)
        {
            row.vec1.push_back(m1.col_ind.at(a1));
            row.vec2.push_back(m1.val.at(a1));
            a1++;
        }
        while (a2 < b2)
        {
            row.vec1.push_back(m2.col_ind.at(a2));
            row.vec2.push_back(m2.val.at(a2));
            a2++;
        }
        v.push_back(row);
        }
        //store the result of the thread in the vector before returning
        result[k] = v;
	};

	//have the threads do the work
	for (size_t i = 0; i < processor_count; ++i){
      	threadPool.push_back(thread(do_work,i));
    }
	//join the threads
	for(auto & val : threadPool){
    	val.join();
	}
    //merge results of each thread into one CSRMatrix
    for (size_t i = 0; i < m1.numRows; i++){
        returnMatrix.col_ind.insert(returnMatrix.col_ind.end(),result[i%processor_count][i/processor_count].vec1.cbegin(),result[i%processor_count][i/processor_count].vec1.cend());
        returnMatrix.val.insert(returnMatrix.val.end(),result[i%processor_count][i/processor_count].vec2.cbegin(),result[i%processor_count][i/processor_count].vec2.cend());
        returnMatrix.row_ptr.push_back(returnMatrix.val.size());
    }
    return returnMatrix;
}

//this matrix code is correct and uses TBB, but it is slow and does not scale well
//it paralyzes over each column
// template<typename T>
// class VectorPair{
//     public:
//     vector<size_t> vec1;
//     vector<T> vec2;
// };
// template<typename T>
// bool compareVectorPairFirstElement(const VectorPair<T>& a, const VectorPair<T>& b) {
//     if (a.vec1.empty() && b.vec1.empty()) {
//         return false;
//     } else if (a.vec1.empty()) {
//         return true;
//     } else if (b.vec1.empty()) {
//         return false;
//     } else {
//         return a.vec1[0] < b.vec1[0];
//     }
// }


// /// @brief Multiplies two compressed sparse row(CSR) matrixes
// /// @exception The number of columns in m1 must equal the number of rows in m2
// /// @tparam T The type of the matrixes
// /// @param m1 The first CSR matrix to multiply
// /// @param m2 The second CSR matrix to multiply
// /// @return The dot product of m1 and m2
// template <typename T>
// CSRMatrix<T> multiply_matrixCSR(CSRMatrix<T> m1, CSRMatrix<T> m2)
// {
//     if (m1.numColumns != m2.numRows)
//     {
//         throw std::invalid_argument("The number of columns in the first matrix must match the number of rows in the second matrix.");
//     }
//     CSRMatrix<T> returnMatrix;
//     returnMatrix.numRows = m1.numRows;
//     returnMatrix.numColumns = m2.numColumns;
//     returnMatrix.row_ptr.push_back(0);
//     CSRMatrix<T> m2t = transpose_matrixCSR(m2);
//     for (size_t i = 0; i < m1.numRows; i++)
//     {
//         tbb::concurrent_vector<VectorPair<T>> columns;
//         tbb::parallel_for( tbb::blocked_range<size_t>(0, m2t.numRows), [&](tbb::blocked_range<size_t> r){
//         VectorPair<T> v;
//         for(size_t j = r.begin(); j < r.end(); j++)
//         {
//             T sum = 0;
//             size_t a1 = m1.row_ptr.at(i);
//             size_t b1 = m1.row_ptr.at(i + 1);
//             size_t a2 = m2t.row_ptr.at(j);
//             size_t b2 = m2t.row_ptr.at(j + 1);
//             while (a1 < b1 && a2 < b2)
//             {
//                 if (m1.col_ind.at(a1) < m2t.col_ind.at(a2))
//                 {
//                     a1++;
//                 }
//                 else if (m1.col_ind.at(a1) > m2t.col_ind.at(a2))
//                 {
//                     a2++;
//                 }
//                 else if (m1.col_ind.at(a1) == m2t.col_ind.at(a2))
//                 {
//                     sum += m1.val.at(a1) * m2t.val.at(a2);
//                     a1++;
//                     a2++;
//                 }
//             }
//             if (sum != 0)
//             {
//                 v.vec1.push_back(j);
//                 v.vec2.push_back(sum);
//             }
//         }
//         columns.push_back(v);
//         });
//         //sort and push on
//         std::sort(columns.begin(),columns.end(),compareVectorPairFirstElement<T>);
//         for(size_t i = 0; i <columns.size();i++){
//             returnMatrix.col_ind.insert(returnMatrix.col_ind.end(),columns[i].vec1.cbegin(),columns[i].vec1.cend());
//             returnMatrix.val.insert(returnMatrix.val.end(),columns[i].vec2.cbegin(),columns[i].vec2.cend());
//         }
//         //cerr << (returnMatrix.val.size()) << endl;
//         returnMatrix.row_ptr.push_back(returnMatrix.val.size());
//     }
//     return returnMatrix;
// }

/// @brief Multiplies two compressed sparse row(CSR) matrixes
/// @exception The number of columns in m1 must equal the number of rows in m2
/// @tparam T The type of the matrixes
/// @param m1 The first CSR matrix to multiply
/// @param m2 The second CSR matrix to multiply
/// @return The dot product of m1 and m2
template <typename T>
CSRMatrix<T> multiply_matrixCSR(CSRMatrix<T> m1, CSRMatrix<T> m2)
{
    if (m1.numColumns != m2.numRows)
    {
        throw std::invalid_argument("The number of columns in the first matrix must match the number of rows in the second matrix.");
    }
    CSRMatrix<T> returnMatrix;
    returnMatrix.numRows = m1.numRows;
    returnMatrix.numColumns = m2.numColumns;
    returnMatrix.row_ptr.push_back(0);
    CSRMatrix<T> m2t = transpose_matrixCSR(m2);
    //set the processor count to the number of available number of threads
    const auto processor_count = std::thread::hardware_concurrency();;
    std::vector<vector<VectorPair<T>>> result;
    //ensures there is a spot for each threads result
    result.resize(processor_count);
    //the vector of threads
  	std::vector<std::thread> threadPool;
	threadPool.reserve(processor_count);
    //lambda to do the multiply for each thread
    auto do_work = [&](int k){
        vector<VectorPair<T>> v;
         //increments by processor count and starts at k to divide work evenly
		for (size_t i = k; i < m1.numRows; i += processor_count){
        VectorPair<T> row;
        for (size_t j = 0; j < m2t.numRows; j++)
        {
            T sum = 0;
            size_t a1 = m1.row_ptr.at(i);
            size_t b1 = m1.row_ptr.at(i + 1);
            size_t a2 = m2t.row_ptr.at(j);
            size_t b2 = m2t.row_ptr.at(j + 1);
            while (a1 < b1 && a2 < b2)
            {
                if (m1.col_ind.at(a1) < m2t.col_ind.at(a2))
                {
                    a1++;
                }
                else if (m1.col_ind.at(a1) > m2t.col_ind.at(a2))
                {
                    a2++;
                }
                else if (m1.col_ind.at(a1) == m2t.col_ind.at(a2))
                {
                    sum += m1.val.at(a1) * m2t.val.at(a2);
                    a1++;
                    a2++;
                }
            }
            if (sum != 0)
            {
                row.vec1.push_back(j);
                row.vec2.push_back(sum);
            }
        }
        v.push_back(row);
        }
        //store the threads result before returning
        result[k] = v;
	};

	//have the threads do the work
	for (size_t i = 0; i < processor_count; ++i){
      	threadPool.push_back(thread(do_work,i));
    }
	//join the threads
	for(auto & val : threadPool){
    	val.join();
	}
    //merge results of each thread into one CSRMatrix
    for (size_t i = 0; i < m1.numRows; i++){
        returnMatrix.col_ind.insert(returnMatrix.col_ind.end(),result[i%processor_count][i/processor_count].vec1.cbegin(),result[i%processor_count][i/processor_count].vec1.cend());
        returnMatrix.val.insert(returnMatrix.val.end(),result[i%processor_count][i/processor_count].vec2.cbegin(),result[i%processor_count][i/processor_count].vec2.cend());
        returnMatrix.row_ptr.push_back(returnMatrix.val.size());
    }
    return returnMatrix;
}

/// @brief gaussian elimination with forward elimination on a square matrix
/// @exception The matrix must be square
/// @param A A dense matrix to eliminate
/// @return returns true if successful
bool gaussian_elimination(std::vector<std::vector<double> > &A) {
    if (A.size() != A[0].size())
    {
        throw std::invalid_argument("The matrix must be square.");
    }
    //Iterate over each row in the matrix
    double pivot;
    timer stopwatch;
    for(size_t i = 0; i < A.size() - 1; i++){
        //Pivot will be the diagonal
        pivot = A[i][i];
        //stopwatch.elapsed();
        //Iterate of the remaining row elements
        tbb::parallel_for( tbb::blocked_range<size_t>(i+1, A[0].size()), [&](tbb::blocked_range<size_t> r){
            for(size_t j = r.begin(); j < r.end(); j++){
                A[i][j] /= pivot;
            }
        });

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
    }
    A[A.size()-1][A[0].size()-1] = 1;

    return true;
}
}

//testing/benchmarking code
// int main() {
//     timer stopwatch;
//     std::vector<vector<double> > parallel;
//     std::vector<vector<double> > serial;
//     parallel.resize(16);
//     serial.resize(16);
//     CSRMatrix<double> m1 = load_fileCSR<double>("stomach.mtx");
//     CSRMatrix<double> m3 = transpose_matrixCSR(m1);
//     // CSRMatrix<double> one = add_matrixCSR(m3,m1);
//     // CSRMatrix<double> two = parallel::add_matrixCSR(m3,m1,16);
//     // if(one.numRows != two.numRows){
//     //     cerr << "yikes" << endl;
//     // }
//     // if(one.numColumns != two.numColumns){
//     //     cerr << "yikes"<< endl;
//     // }
//     // if(one.row_ptr.size() != two.row_ptr.size()){
//     //     cerr << "yikes"<< endl;
//     // }
//     // if(one.col_ind.size() != two.col_ind.size()){
//     //     cerr << "yikes"<< endl;
//     // }
//     // if(one.val.size() != two.val.size()){
//     //     cerr << "yikes"<< endl;
//     // }
//     // for(size_t i = 0 ; i <one.row_ptr.size();i++){
//     //     if(one.row_ptr[i] != two.row_ptr[i]){
//     //     cerr << "yikes"<< endl;
//     // }
//     // }
//     // for(size_t i = 0 ; i <one.col_ind.size();i++){
//     //      if(one.col_ind[i] != two.col_ind[i]){
//     //     cerr << "yikes"<< endl;
//     // }
//     // }
//     // for(size_t i = 0 ; i <one.val.size();i++){
//     //     if(one.val[i] != two.val[i]){
//     //     cerr << "yikes"<< endl;
//     // }
//     // }
//     // cerr<< "yay!";
//     for(size_t i = 0 ; i <16;i++){
//         for(size_t j = 0; j <10;j++){
//             stopwatch.elapsed();
//             parallel::add_matrixCSR(m1,m3,i+1);
//             parallel[i].push_back(stopwatch.elapsed());
//             stopwatch.elapsed();
//             add_matrixCSR(m1,m3);
//             serial[i].push_back(stopwatch.elapsed());
//         }
//     }
//     for(size_t i = 0; i < parallel.size();i++){
//         double time = 0;
//         for(size_t j = 0 ; j < parallel[0].size();j++){
//             time += parallel[i][j];
//         }
//         cerr<< time/parallel[0].size() << ",";
//     }
//     cerr<< endl;
//     for(size_t i = 0; i < serial.size();i++){
//         double time = 0;
//         for(size_t j = 0 ; j < serial[0].size();j++){
//             time += serial[i][j];
//         }
//         cerr<< time/serial[0].size() << ",";
//     }
//     cerr<< endl;
//     return 0;
// }