#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

template<typename T>
	class CSRMatrix{
        // struct colValStruct {
        //     int col;
        //     T value;
        // };
        public:
        size_t numRows, numColumns;
        // vector<colValStruct*> rowPointersVector;
        //fist value is col_ind and the second value is the value
        // vector<colValStruct> columnValueVector;
        vector<T> val;
        vector<size_t> col_ind;
        vector<size_t> row_ptr;
    };
    //loadfile
    //savefile
template<typename T>
    T get_matrixCSR(CSRMatrix<T> m1, size_t row, size_t col){
        if(m1.numRows <= row){
            //throw error
        }
        if(m1.numColumns <= col){
            //throw error
        }
        for (size_t i = m1.row_ptr.at(row); i < m1.row_ptr.at(row+1); i++){
            if(m1.col_ind.at(i) == col){
                return m1.val.at(i);
            }
        }
        return 0;
    }

template<typename T>
    CSRMatrix<T> from_vector(vector<vector<T> > &array){
        CSRMatrix<T> returnMatrix;
        returnMatrix.numRows = array.size();
        returnMatrix.numColumns = array.at(0).size();
        returnMatrix.row_ptr.push_back(0);
        for (size_t i = 0; i < array.size(); i++){
            for (size_t j = 0; j < array.at(i).size(); j++){
                if(array.at(i).at(j) != 0){
                    returnMatrix.val.push_back(array.at(i).at(j));
                    returnMatrix.col_ind.push_back(j);
                }
            }
            returnMatrix.row_ptr.push_back(returnMatrix.val.size());
        }
        return returnMatrix;
    }

template<typename T>
    void print_matrixCSR(CSRMatrix<T> m1){
        for (size_t i = 0; i < m1.numRows; i++) {
            for (size_t j = 0; j < m1.numColumns; j++) {
                cout << get_matrixCSR(m1, i, j) << " ";
            }
            cout << endl;
        }
    }

template<typename T>
    CSRMatrix<T> add_matrixCSR(CSRMatrix<T> m1, CSRMatrix<T> m2){
        if(m1.numRows!= m2.numRows){
            throw std::invalid_argument("The number of rows in the first matrix must match the number of rows in the second matrix.");
        }
        if(m1.numColumns!= m2.numColumns){
            throw std::invalid_argument("The number of columns in the first matrix must match the number of columns in the second matrix.");
        }
        CSRMatrix<T> returnMatrix;
        returnMatrix.numRows = m1.numRows;
        returnMatrix.numColumns = m1.numColumns;
        returnMatrix.row_ptr.push_back(0);
        
        for (size_t i = 0; i < m1.numRows; i++) {
            size_t a1 = m1.row_ptr.at(i);
            size_t b1 = m1.row_ptr.at(i+1);
            size_t a2 = m2.row_ptr.at(i);
            size_t b2 = m2.row_ptr.at(i+1);
            while(a1 < b1 && a2 < b2){
                if (m1.col_ind.at(a1) < m2.col_ind.at(a2)) {
                    returnMatrix.val.push_back(m1.val.at(a1));
                    returnMatrix.col_ind.push_back(m1.col_ind.at(a1));
                    a1++;
                } else if (m1.col_ind.at(a1) > m2.col_ind.at(a2)) {
                    returnMatrix.val.push_back(m2.val.at(a2));
                    returnMatrix.col_ind.push_back(m2.col_ind.at(a2));
                    a2++;
                }
                else if (m1.col_ind.at(a1) == m2.col_ind.at(a2)) {
                    T value = m1.val.at(a1) + m2.val.at(a2);
                    if (value != 0) {
                        returnMatrix.val.push_back(value);
                        returnMatrix.col_ind.push_back(m1.col_ind.at(a1));
                    }
                    a1++;
                    a2++;
                }
            }
            while (a1 < b1) {
                returnMatrix.val.push_back(m1.val.at(a1));
                returnMatrix.col_ind.push_back(m1.col_ind.at(a1));
                a1++;
            }
            while (a2 < b2) {
                returnMatrix.val.push_back(m2.val.at(a2));
                returnMatrix.col_ind.push_back(m2.col_ind.at(a2));
                a2++;
            }
            returnMatrix.row_ptr.push_back(returnMatrix.val.size());
        }
        return returnMatrix;
    }

template<typename T>
    CSRMatrix<T> transpose_matrixCSR(CSRMatrix<T> m1){
        CSRMatrix<T> returnMatrix;
        returnMatrix.numRows = m1.numColumns;
        returnMatrix.numColumns = m1.numRows;
        returnMatrix.row_ptr.push_back(0);
        vector<size_t> row_count(m1.numColumns, 0);
        returnMatrix.val = vector<T>(m1.val.size());
        returnMatrix.col_ind = vector<size_t>(m1.col_ind.size());
        for (size_t i = 0; i < m1.numRows; i++) {
            for (size_t j = m1.row_ptr.at(i); j < m1.row_ptr.at(i+1); j++) {
                row_count.at(m1.col_ind.at(j))++;
            }
        }
        for (size_t i = 0; i < m1.numColumns; i++) {
            returnMatrix.row_ptr.push_back(returnMatrix.row_ptr.at(i)+row_count.at(i));
        }
        row_count = vector<size_t>(m1.numColumns, 0);
        for (size_t i = 0; i < m1.numRows; i++) {
            for (size_t j = m1.row_ptr.at(i); j < m1.row_ptr.at(i+1); j++) {
                size_t index = returnMatrix.row_ptr.at(m1.col_ind.at(j)) + row_count.at(m1.col_ind.at(j));
                returnMatrix.val.at(index) = m1.val.at(j);
                returnMatrix.col_ind.at(index) = i;
                row_count.at(m1.col_ind.at(j))++;
            }
        }
        return returnMatrix;
    }

template<typename T>
    CSRMatrix<T> multiply_matrixCSR(CSRMatrix<T> m1, CSRMatrix<T> m2) {
        if(m1.numColumns!= m2.numRows){
            throw std::invalid_argument("The number of columns in the first matrix must match the number of rows in the second matrix.");
        }
        CSRMatrix<T> returnMatrix;
        returnMatrix.numRows = m1.numRows;
        returnMatrix.numColumns = m2.numColumns;
        returnMatrix.row_ptr.push_back(0);
        CSRMatrix<T> m2t = transpose_matrixCSR(m2);
        for (size_t i = 0; i < m1.numRows; i++) {
            for (size_t j = 0; j < m2t.numRows; j++) {
                T sum = 0;
                size_t a1 = m1.row_ptr.at(i);
                size_t b1 = m1.row_ptr.at(i+1);
                size_t a2 = m2t.row_ptr.at(j);
                size_t b2 = m2t.row_ptr.at(j+1);
                while(a1 < b1 && a2 < b2){
                    if (m1.col_ind.at(a1) < m2t.col_ind.at(a2)) {
                        a1++;
                    } else if (m1.col_ind.at(a1) > m2t.col_ind.at(a2)) {
                        a2++;
                    }
                    else if (m1.col_ind.at(a1) == m2t.col_ind.at(a2)) {
                        sum += m1.val.at(a1)*m2t.val.at(a2);
                        a1++;
                        a2++;
                    }
                }
                if (sum != 0) {
                    returnMatrix.val.push_back(sum);
                    returnMatrix.col_ind.push_back(j);
                }
            }
            returnMatrix.row_ptr.push_back(returnMatrix.val.size());
        }
        return returnMatrix;
    }

// int main() {
//     vector<vector<int> > array = vector<vector<int> >(3, vector<int>(3));
//     array[0][0] = 1;
//     array[0][1] = 2;
//     array[0][2] = 3;
//     array[1][0] = 4;
//     array[1][1] = 5;
//     array[1][2] = 6;
//     array[2][0] = 7;
//     array[2][1] = 8;
//     array[2][2] = 9;
//     // array[3][0] = 10;
//     // array[3][1] = 11;
//     // array[3][2] = 12;

//     vector<vector<int> > array2 = vector<vector<int> >(3, vector<int>(3));
//     array2[0][0] = 9;
//     array2[0][1] = -4;
//     array2[0][2] = 7;
//     array2[1][0] = -6;
//     array2[1][1] = -5;
//     array2[1][2] = 4;
//     array2[2][0] = 3;
//     array2[2][1] = 2;
//     array2[2][2] = 1;
//     // array2[3][0] = 10;
//     // array2[3][1] = 11;
//     // array2[3][2] = 12;

//     CSRMatrix<int> m1 = from_vector(array);
//     CSRMatrix<int> m2 = from_vector(array2);
//     CSRMatrix<int> m3 = add_matrixCSR(m1, transpose_matrixCSR(m2));
//     print_matrixCSR(m3);
//     m3 = transpose_matrixCSR(m3);
//     // m3 = transpose_matrixCSR(m3);
//     print_matrixCSR(m3);

//     return 0;
// }