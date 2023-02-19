#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

template<typename T>
	class CSRMatrix{
        public:
        int numRows, numColumns;
        vector<pair*> rowPointersVector;
        //fist value is col_ind and the second value is the value
        vector<pair<int,T>> columnValueVector;
    };
    //loadfile
    //savefile

    double get_matrixCSR(CSRMatrix<double> m1, int row, int col){
        if(m1.numRows< row){
            //throw error
        }
        if(m1.numColumns < col){
            //throw error
        }
        //if the value being looked for is in the last row
        if(m1.numRows == row-1){
            for(auto it = m1.rowPointersVector.at(row); it != m1.columnValueVector.end();it++){
                if(it.first == col){
                    return it.second;
                }
                return -1;
            }
        }
        for(auto it = m1.rowPointersVector.at(row); it != m1.rowPointersVector.at(row+1);it++){
            if(it.first == col){
                return it.second;
            }
        }
        return -1;

    }

    CSRMatrix<double> add_matrixCSR(CSRMatrix<double> m1, CSRMatrix<double> m2){
        if(m1.numRows!= m2.numRows){
            //throw error
        }
        if(m1.numColumns!= m2.numColumns){
            //throw error
        }
        CSRMatrix<double> returnMatrix;
        for(size_t i =0; i < m2.rowPointersVector.size();i++){

        }

    }
    //TODO: Multiply and transpose