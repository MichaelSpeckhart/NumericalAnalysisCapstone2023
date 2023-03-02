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
        // for(auto it = m1.rowPointersVector.at(row); it != m1.rowPointersVector.at(row+1);it++){
        //     if(it->col == col){
        //         return it->value;
        //     }
        // }
        // return -1;
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
    CSRMatrix<T> add_matrixCSR(CSRMatrix<T> m1, CSRMatrix<T> m2){
        if(m1.numRows!= m2.numRows){
            //throw error
        }
        if(m1.numColumns!= m2.numColumns){
            //throw error
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
                    returnMatrix.val.push_back(m1.val.at(a1)+m2.val.at(a2));
                    returnMatrix.col_ind.push_back(m1.col_ind.at(a1));
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

    //TODO: Multiply and transpose

int main() {
    vector<vector<int> > array = vector<vector<int> >(3, vector<int>(3));
    array[0][0] = 1;
    // array[0][1] = 2;
    // array[0][2] = 3;
    array[1][0] = 4;
    array[1][1] = 5;
    array[1][2] = 6;
    // array[2][0] = 7;
    array[2][1] = 8;
    array[2][2] = 9;

    vector<vector<int> > array2 = vector<vector<int> >(3, vector<int>(3));
    // array2[0][0] = 1;
    array2[0][1] = 2;
    array2[0][2] = 3;
    // array2[1][0] = 4;
    array2[1][1] = -5;
    // array2[1][2] = 6;
    array2[2][0] = 7;
    array2[2][1] = 8;
    // array2[2][2] = 9;

    CSRMatrix<int> m1 = from_vector<int>(array);
    CSRMatrix<int> m2 = from_vector<int>(array2);
    CSRMatrix<int> m3 = add_matrixCSR<int>(m1, m2);
    cout << get_matrixCSR(m3, 0, 0) << endl;
    cout << get_matrixCSR(m3, 0, 1) << endl;
    cout << get_matrixCSR(m3, 0, 2) << endl;
    cout << get_matrixCSR(m3, 1, 0) << endl;
    cout << get_matrixCSR(m3, 1, 1) << endl;
    cout << get_matrixCSR(m3, 1, 2) << endl;
    cout << get_matrixCSR(m3, 2, 0) << endl;
    cout << get_matrixCSR(m3, 2, 1) << endl;
    cout << get_matrixCSR(m3, 2, 2) << endl;
    return 0;
}